#include <typeinfo>
#include "grid_node.h"
#include "indexer_2d.h"
#include "indexer_2d_by2.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int GridNode::max_refine_level = 5;
Real GridNode::min_refine_density = 1.0e-6;
Real GridNode::refine_alpha = 1.0;

Real GridNode::t_norm = 1.0;
void GridNode::set_time_normal(Real t) {
	t_norm = t;

}

void GridNode::set_origin(const _3Vec& O) {
	Grid::set_origin(O);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->set_origin(O);
		}
	}
}

void GridNode::write(FILE* fp) const {
	bool c[OCT_NCHILD];
	Grid::write(fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		c[i] = (get_child(i) != NULL);
	}
	fwrite(&c, sizeof(bool), 8, fp);
	fwrite(&refine_alpha, sizeof(refine_alpha), 1, fp);
	fwrite(&age, sizeof(age), 1, fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (c[i]) {
			static_cast<const GridNode*> (get_child(i))->write(fp);
		}
	}
}

void GridNode::read(FILE* fp) {
	bool c[OCT_NCHILD];
	Grid::read(fp);
	fread(c, sizeof(bool), 8, fp);
	fread(&refine_alpha, sizeof(refine_alpha), 1, fp);
	fread(&age, sizeof(age), 1, fp);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (c[i]) {
			this->create_child(i);
			static_cast<GridNode*> (get_child(i))->read(fp);
		}
	}
}

void GridNode::enforce_boundaries() {
	const Indexer2d indexer(bw, GNX - 1 - bw, bw, GNX - 1 - bw);
	int i, j, k;
	const GridInterface* sib;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i,j,k,sib)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		sib = dynamic_cast<const GridInterface*> (get_sibling(XL));
		assert( sib != NULL);
		for (i = 0; i < bw; i++) {
			(*this)(i, j, k) = (*sib)(GNX - 2 * bw + i, j, k);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(XU));
		assert( sib != NULL);
		for (i = GNX - bw; i < GNX; i++) {
			(*this)(i, j, k) = (*sib)(-GNX + 2 * bw + i, j, k);
		}
		i = indexer.x(index);
		k = indexer.y(index);
		sib = dynamic_cast<const GridInterface*> (get_sibling(YL));
		assert( sib != NULL);
		for (j = 0; j < bw; j++) {
			(*this)(i, j, k) = (*sib)(i, GNX - 2 * bw + j, k);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(YU));
		assert( sib != NULL);
		for (j = GNX - bw; j < GNX; j++) {
			(*this)(i, j, k) = (*sib)(i, -GNX + 2 * bw + j, k);
		}
		i = indexer.x(index);
		j = indexer.y(index);
		sib = dynamic_cast<const GridInterface*> (get_sibling(ZL));
		assert( sib != NULL);
		for (k = 0; k < bw; k++) {
			(*this)(i, j, k) = (*sib)(i, j, GNX - 2 * bw + k);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(ZU));
		assert( sib != NULL);
		for (k = GNX - bw; k < GNX; k++) {
			(*this)(i, j, k) = (*sib)(i, j, -GNX + 2 * bw + k);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->enforce_boundaries();
		}
	}
#ifdef RECON3D
	if( get_level() == 0 ) {
		enforce_edge_boundaries();
	}
#endif
}

void GridNode::enforce_edge_boundaries() {
	int i, j, k;
	const GridInterface* sib;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i,k,sib)
	for (j = bw; j < GNX - bw; j++) {
		sib = dynamic_cast<const GridInterface*> (get_sibling(XL));
		assert( sib != NULL);
		for (i = 0; i < bw; i++) {
			k = 0;
			(*this)(i, j, k) = (*sib)(GNX - 2 * bw + i, j, k);
			(*this)(i, k, j) = (*sib)(GNX - 2 * bw + i, k, j);
			k = GNX - BW - 1;
			(*this)(i, j, k) = (*sib)(GNX - 2 * bw + i, j, k);
			(*this)(i, k, j) = (*sib)(GNX - 2 * bw + i, k, j);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(XU));
		assert( sib != NULL);
		for (i = GNX - bw; i < GNX; i++) {
			k = 0;
			(*this)(i, j, k) = (*sib)(-GNX + 2 * bw + i, j, k);
			(*this)(i, k, j) = (*sib)(-GNX + 2 * bw + i, k, j);
			k = GNX - BW - 1;
			(*this)(i, j, k) = (*sib)(-GNX + 2 * bw + i, j, k);
			(*this)(i, k, j) = (*sib)(-GNX + 2 * bw + i, k, j);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(YL));
		assert( sib != NULL);
		for (i = 0; i < bw; i++) {
			k = 0;
			(*this)(j, i, k) = (*sib)(j, GNX - 2 * bw + i, k);
			k = GNX - BW - 1;
			(*this)(j, i, k) = (*sib)(j, GNX - 2 * bw + i, k);
		}
		sib = dynamic_cast<const GridInterface*> (get_sibling(YU));
		assert( sib != NULL);
		for (i = GNX - bw; i < GNX; i++) {
			k = 0;
			(*this)(j, i, k) = (*sib)(j, -GNX + 2 * bw + i, k);
			k = GNX - BW - 1;
			(*this)(j, i, k) = (*sib)(j, -GNX + 2 * bw + i, k);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->enforce_boundaries();
		}
	}
}

Real GridNode::max_dt() const {
	Real dtinv = 0.0;
	Real cdtinv;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			cdtinv = 1.0 / static_cast<const GridNode*> (get_child(i))->max_dt();
			dtinv = max(cdtinv, dtinv);
		}
	}
	dtinv = max(1.0 / Grid::max_dt(), dtinv);
	return 1.0 / dtinv;
}

void GridNode::debug() const {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			if ((*this)(i, j, k).rho() < 0.0) {
				printf("%e %i %i %i %i", (*this)(i, j, k).rho(), i, j, k, get_level());
				printf("%i %i %i", get_location()[0], get_location()[1], get_location()[2]);
				abort();
			}
		}
	}
}

void GridNode::inject_from_parent(ChildIndex c) {
	const GridNode* p = static_cast<const GridNode*> (get_parent());
	const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	Vector<Real, STATE_NF> s1, s2, s3;
	State u;
	int k, j, k0, j0, i, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(s1,s2,s3,u,k,j,i,k0,j0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		k0 = (bw + k) / 2 + c.get_z() * (GNX / 2 - bw);
		j0 = (bw + j) / 2 + c.get_y() * (GNX / 2 - bw);
		for (i = bw, i0 = bw + c.get_x() * (GNX / 2 - bw); i < GNX - bw; i += 2, i0++) {
			u = (*p)(i0, j0, k0);
#ifdef PIECEWISE_CONSTANT_INJECTION
			(*this)(i + 0, j + 0, k + 0) = (*this)(i + 1, j + 0, k + 0) = (*this)(i + 0, j + 1, k + 0) = (*this)(i + 1,
					j + 1, k + 0) = (*this)(i + 0, j + 0, k + 1) = (*this)(i + 1, j + 0, k + 1) = (*this)(i + 0, j + 1,
					k + 1) = (*this)(i + 1, j + 1, k + 1) = u;
#else
			s1 = minmod((*p)(i0 + 1, j0, k0) - u, u - (*p)(i0 - 1, j0, k0));
			s2 = minmod((*p)(i0, j0 + 1, k0) - u, u - (*p)(i0, j0 - 1, k0));
			s3 = minmod((*p)(i0, j0, k0 + 1) - u, u - (*p)(i0, j0, k0 - 1));
			(*this)(i + 0, j + 0, k + 0) = u - (s1 + s2 + s3) * 0.25;
			(*this)(i + 1, j + 0, k + 0) = u + (s1 - s2 - s3) * 0.25;
			(*this)(i + 0, j + 1, k + 0) = u - (s1 - s2 + s3) * 0.25;
			(*this)(i + 1, j + 1, k + 0) = u + (s1 + s2 - s3) * 0.25;
			(*this)(i + 0, j + 0, k + 1) = u - (s1 + s2 - s3) * 0.25;
			(*this)(i + 1, j + 0, k + 1) = u + (s1 - s2 + s3) * 0.25;
			(*this)(i + 0, j + 1, k + 1) = u - (s1 - s2 - s3) * 0.25;
			(*this)(i + 1, j + 1, k + 1) = u + (s1 + s2 + s3) * 0.25;
#endif
		}
	}
}

void GridNode::error_from_parent() {
#ifdef DEN_REF
	return;
#endif
	if (get_level() != 0) {
		const ChildIndex c = my_child_index();
		const GridNode* p = static_cast<const GridNode*> (get_parent());
		const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
		State u;
		int k, j, k0, j0, i, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(u,k,j,i,k0,j0,i0)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			k0 = (bw + k) / 2 + c.get_z() * (GNX / 2 - bw);
			j0 = (bw + j) / 2 + c.get_y() * (GNX / 2 - bw);
			for (i = bw, i0 = bw + c.get_x() * (GNX / 2 - bw); i < GNX - bw; i += 2, i0++) {
				u = p->differential(i0, j0, k0);
				E(i + 0, j + 0, k + 0) = u - differential(i + 0, j + 0, k + 0);
				E(i + 1, j + 0, k + 0) = u - differential(i + 1, j + 0, k + 0);
				E(i + 0, j + 1, k + 0) = u - differential(i + 0, j + 1, k + 0);
				E(i + 1, j + 1, k + 0) = u - differential(i + 1, j + 1, k + 0);
				E(i + 0, j + 0, k + 1) = u - differential(i + 0, j + 0, k + 1);
				E(i + 1, j + 0, k + 1) = u - differential(i + 1, j + 0, k + 1);
				E(i + 0, j + 1, k + 1) = u - differential(i + 0, j + 1, k + 1);
				E(i + 1, j + 1, k + 1) = u - differential(i + 1, j + 1, k + 1);
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->GridNode::error_from_parent();
		}
	}
}

void GridNode::create_child(const ChildIndex& c) {
	GridNode* ptr;
	OctNode* sib;
	OctNode::create_child(c);
	Vector<int, 3> child_offset;
	child_offset = (get_offset() * 2 + bw) + c.vector() * (GNX - 2 * bw);
	ptr = static_cast<GridNode*> (get_child(c));
	ptr->set_time(get_time());
	ptr->set_dx(get_dx() * 0.5);
	ptr->set_offset(child_offset);
	ptr->set_origin(get_origin());
	for (int i = 0; i < OCT_NSIB; i++) {
		if (ptr->get_sibling(i) == NULL) {
			if (ptr->is_phys_bound(OctFace(i))) {
				sib = this->new_phys_bound();
			} else {
				sib = this->new_amr_bound();
			}
			ptr->set_sibling(OctFace(i), sib);
			sib->set_sibling(invert(OctFace(i)), ptr);
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		ptr->refine[i] = false;
	}
	ptr->inject_from_parent(c);
	//	ptr->enforce_boundaries();
	assert( (ptr->get_sibling(XL)!= NULL));
	assert( (ptr->get_sibling(YL)!= NULL));
	assert( (ptr->get_sibling(ZL)!= NULL));
	assert( (ptr->get_sibling(XU)!= NULL));
	assert( (ptr->get_sibling(YU)!= NULL));
	assert( (ptr->get_sibling(ZU)!= NULL));
}

GridNode* GridNode::new_octnode() const {
	return new GridNode;
}

GridAMRBound* GridNode::new_amr_bound() const {
	return new GridAMRBound();
}

GridPhysBound* GridNode::new_phys_bound() const {
	return new GridPhysBound();
}

void GridNode::destroy_child(const ChildIndex& c) {
	OctNode* sib[OCT_NSIB];
	OctNode* amr;
	for (int i = 0; i < OCT_NSIB; i++) {
		sib[i] = get_child(c)->get_sibling(i);
	}
	OctNode::destroy_child(c);
	for (int i = 0; i < OCT_NSIB; i++) {
		assert (sib[i] != NULL);
		if (sib[i]->is_real()) {
			amr = this->new_amr_bound();
			sib[i]->set_sibling(invert(OctFace(i)), amr);
			amr->set_sibling(OctFace(i), sib[i]);
		} else {
			delete sib[i];
		}
	}
}

void GridNode::add_difs(Real dt, Real beta) {
	age++;
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
#ifndef DEN_REF
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			E0(i, j, k) = (E0(i, j, k) + E(i, j, k)) * beta;
		}
	}
#endif
	Grid::add_difs(dt, beta);
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->GridNode::add_difs(dt, beta);
		}
	}
}

void GridNode::clear_difs() {
	Grid::clear_difs();
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->clear_difs();
		}
	}

}

void GridNode::compute_x_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->compute_x_flux();
		}
	}
	Grid::compute_x_flux();
}

void GridNode::compute_y_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->compute_y_flux();
		}
	}
	Grid::compute_y_flux();
}

void GridNode::compute_z_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->compute_z_flux();
		}
	}
	Grid::compute_z_flux();
}

void GridNode::restore() {
	Grid::restore();
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->restore();
		}
	}
}

void GridNode::store() {
	int i;
#ifndef DEN_REF
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	int k, j;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			E0(i, j, k) = Vector<Real, STATE_NF> (0.0);
		}
	}
#endif
	Grid::store();
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->store();
		}
	}
}

void GridNode::sum_x_difs() {
	Grid::sum_x_difs();
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->sum_x_difs();
		}
	}
	if (get_level() == 0) {
		sum_x_flow_off();
	}
}

void GridNode::sum_y_difs() {
	Grid::sum_y_difs();
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->sum_y_difs();
		}
	}
	if (get_level() == 0) {
		sum_y_flow_off();
	}
}

void GridNode::sum_z_difs() {
	Grid::sum_z_difs();
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->sum_z_difs();
		}
	}
	if (get_level() == 0) {
		sum_z_flow_off();
	}
}

void GridNode::set_time(Real t) {
	if (t == 0.0) {
		age = -1;
	}
	Grid::set_time(t);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->set_time(t);
		}
	}
}

void GridNode::set_dx(Real a) {
	Grid::set_dx(a);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->set_dx(a);
		}
	}
}

int GridNode::refine_flag_cnt() const {
	int cnt = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (refine[i]) {
			cnt++;
		}
		if (this->get_child(i) != NULL) {
			cnt += static_cast<const GridNode*> (get_child(i))->refine_flag_cnt();
		}
	}
	return cnt;
}

int GridNode::fine_point_count() const {
	if (get_level() <= GridNode::max_refine_level) {
		int cnt = 0;
		const int cnt_per = ((GNX - 2 * BW) * (GNX - 2 * BW) * (GNX - 2 * BW)) / 8;
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) == NULL) {
				if (refine[i]) {
					cnt += 8 * cnt_per;
				} else {
					cnt += cnt_per;
				}
			}
			if (this->get_child(i) != NULL) {
				cnt += static_cast<const GridNode*> (get_child(i))->fine_point_count();
			}
		}
		return cnt;
	} else {
		return 0;
	}
}

void GridNode::exec_function(void(*f)(GridNode*), int minlev) {
	if (get_level() >= minlev) {
		f(this);
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->exec_function(f, minlev);
		}
	}
}

void GridNode::clear_refine_flags() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		refine[i] = false;
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->clear_refine_flags();
		}
	}
}

void GridNode::set_refine_flags() {
	if (get_level() < MINLEVEL) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			refine[i] = true;
		}
	} else if (get_level() < GridNode::max_refine_level) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				const GridNode* child = static_cast<GridNode*> (get_child(i));
				for (int j = 0; j < OCT_NCHILD; j++) {
					if (child->refine[j] || (child->get_child(j) != NULL)) {
						refine[i] = true;
					}
				}
			}
		}
		const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
		ChildIndex c;
		Real refine_criteria, e, tot;
		int k, j, i;
		bool test;
#pragma omp parallel for private(c,e,k,j,i,test,refine_criteria,tot) schedule(dynamic,1)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = bw; i < GNX - bw; i++) {
				c.set_x(2 * i / GNX);
				c.set_y(2 * j / GNX);
				c.set_z(2 * k / GNX);
				bool test;
				test = (*this)(i, j, k).rho() > GridNode::min_refine_density;
#pragma omp critical
				{
					refine[c] = refine[c] || test;
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			if (get_child(i)->nchildren() > 0) {
				refine[i] = true;
			}
			static_cast<GridNode*> (get_child(i))->set_refine_flags();
			const int age_c = static_cast<const GridNode*> (get_child(i))->age;
			if (age_c < MIN_GRID_LIFE && get_time() != 0.0 && get_level() < GridNode::max_refine_level) {
				refine[i] = true;
			}
		}
	}
}

void GridNode::enforce_proper_nesting() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL && this->refine[i]) {
			static_cast<GridNode*> (get_child(i))->enforce_proper_nesting();
		}
	}
	if (get_level() > 0) {
		const ChildIndex mi = my_child_index();
		OctNode* parent = get_parent();
		ChildIndex ci, i;
		GridNode* g;
		OctNode* tmp;
#pragma omp parallel for private(ci,i,g,tmp) schedule(dynamic,1)
		for (int ci0 = 0; ci0 < OCT_NCHILD; ci0++) {
			ci = ci0;
			if (refine[ci]) {
				i = mi;
				i.flip_x();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					g = dynamic_cast<GridNode*> (parent->get_sibling(ci.x_face()));
					assert(g!=NULL);
				} else {
					g = dynamic_cast<GridNode*> (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_y();
				g = NULL;
				if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
					g = dynamic_cast<GridNode*> (parent->get_sibling(ci.y_face()));
					assert(g!=NULL);
				} else {
					g = dynamic_cast<GridNode*> (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_z();
				g = NULL;
				if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
					g = dynamic_cast<GridNode*> (parent->get_sibling(ci.z_face()));
					assert(g!=NULL);
				} else {
					g = dynamic_cast<GridNode*> (parent);
					assert(g!=NULL);
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_y();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						tmp = parent->get_sibling(ci.x_face());
						tmp = tmp->get_sibling(ci.y_face());
						g = dynamic_cast<GridNode*> (tmp);
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.x_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.y_face()));
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_z();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						tmp = parent->get_sibling(ci.x_face());
						tmp = tmp->get_sibling(ci.z_face());
						g = dynamic_cast<GridNode*> (tmp);
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.x_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.z_face()));
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_y();
				i.flip_z();
				g = NULL;
				if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						tmp = parent->get_sibling(ci.y_face());
						tmp = tmp->get_sibling(ci.z_face());
						g = dynamic_cast<GridNode*> (tmp);
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.y_face()));
						assert(g!=NULL);
					}
				} else {
					if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
						g = dynamic_cast<GridNode*> (parent->get_sibling(ci.z_face()));
						assert(g!=NULL);
					} else {
						g = dynamic_cast<GridNode*> (parent);
						assert(g!=NULL);
					}
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
				i = mi;
				i.flip_x();
				i.flip_y();
				i.flip_z();
				g = NULL;
				if ((mi.get_x() == ci.get_x()) && !is_phys_bound(ci.x_face())) {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.x_face());
							assert(tmp!=NULL);
							tmp = tmp->get_sibling(ci.y_face());
							assert(tmp!=NULL);
							tmp = tmp->get_sibling(ci.z_face());
							assert(tmp!=NULL);
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
							if (g == NULL) {
								printf("%s\n", typeid(*tmp).name());
							}
						} else {
							tmp = parent->get_sibling(ci.x_face());
							tmp = tmp->get_sibling(ci.y_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						}
					} else {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.x_face());
							tmp = tmp->get_sibling(ci.z_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						} else {
							tmp = parent->get_sibling(ci.x_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						}
					}
				} else {
					if ((mi.get_y() == ci.get_y()) && !is_phys_bound(ci.y_face())) {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.y_face());
							tmp = tmp->get_sibling(ci.z_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						} else {
							tmp = parent->get_sibling(ci.y_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						}
					} else {
						if ((mi.get_z() == ci.get_z()) && !is_phys_bound(ci.z_face())) {
							tmp = parent->get_sibling(ci.z_face());
							g = dynamic_cast<GridNode*> (tmp);
							assert(g!=NULL);
						} else {
							g = dynamic_cast<GridNode*> (parent);
							assert(g!=NULL);
						}
					}
				}
				if (g != NULL) {
					assert( g->is_real());
#pragma omp critical
					g->refine[i] = true;
				}
			}
		}
	}
}

void GridNode::propagate_refine_flags_up() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->propagate_refine_flags_up();
			for (int j = 0; j < OCT_NCHILD; j++) {
				if (static_cast<GridNode*> (get_child(i))->refine[j]) {
					refine[i] = true;
					break;
				}
			}
		}
	}
}

bool GridNode::use_refine_flags() {
	bool rc = false;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (refine[i] && (get_child(i) == NULL) && (get_level() + 1 <= max_refine_level)) {
			create_child(i);
			rc = true;
		} else if (!refine[i] && (get_child(i) != NULL)) {
			destroy_child(i);
			rc = true;
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			bool test = static_cast<GridNode*> (get_child(i))->use_refine_flags();
			rc = rc || test;
		}
	}
	return rc;
}

void GridNode::mirror_refine(const GridNode* mirror, int dir) {
	if (get_level() < GridNode::max_refine_level) {
		int ir;
		const int mask = 1 << dir;
		const GridNode* child_mirror;
		for (int i = 0; i < OCT_NCHILD; i++) {
			ir = i ^ mask;
			if (mirror->refine[ir]) {
				refine[i] = true;
			}
		}
		for (int i = 0; i < OCT_NCHILD; i++) {
			ir = i ^ mask;
			child_mirror = static_cast<const GridNode*> (mirror->get_child(ir));
			if (get_child(i) != NULL && child_mirror != NULL) {
				static_cast<GridNode*> (get_child(i))->mirror_refine(child_mirror, dir);
			}
		}
	}
}

bool GridNode::check_for_refine(bool clear) {
	bool rc = false;
	int inc, dec;
	if (max_refine_level == 0) {
		return false;
	}
	//	int dir = 0;
	//	int iters = 0;
	bool done = false;
	inc = dec = 0;
	while (!done) {
		int npoints;
		if (clear) {
			clear_refine_flags();
		}
		set_refine_flags();
		propagate_refine_flags_up();
		enforce_proper_nesting();
#ifdef MIRROR_REFINE_X
		mirror_refine(this, 0);
#endif
#ifdef MIRROR_REFINE_Y
		mirror_refine(this, 1);
#endif
#ifdef MIRROR_REFINE_Z
		mirror_refine(this,2);
#endif
		propagate_refine_flags_up();
		npoints = this->fine_point_count();
		/*	if (npoints > MAX_POINTS && inc < 1) {
		 printf("Fine point count of %i exceeds limit of %i\n", npoints, MAX_POINTS);
		 GridNode::max_refine_level--;
		 printf("Decreasing maximum refinement level to %i\n", GridNode::max_refine_level);
		 inc++;
		 } else if (dec < 1 && npoints < MAX_POINTS / 2 && this->max_level() >= GridNode::max_refine_level) {
		 printf("Fine point count of %i violates lower limit of %i\n", npoints, MAX_POINTS / 2);
		 GridNode::max_refine_level++;
		 printf("Increasing maximum refinement level to %i\n", GridNode::max_refine_level);
		 dec++;
		 } else {
		 done = true;

		 }*/
		break;
	}
	/*	do {
	 iters++;
	 clear_refine_flags();
	 set_refine_flags();
	 propagate_refine_flags_up();
	 enforce_proper_nesting();
	 propagate_refine_flags_up();
	 int maxlgrids = ((1 << (6 * max_level() + 6)) - 1) / 7;
	 int ng = refine_flag_cnt() + 1;
	 //		printf("%i %i\n", ng, maxlgrids);
	 if (iters > 2048) {
	 done = true;
	 } else if (ng > MAX_GRIDS && dir != 1) {
	 dir = -1;
	 refine_alpha *= 1.1;
	 done = false;
	 //			printf("Decreasing refinement %i %e\n", refine_flag_cnt() + 1, refine_alpha);
	 } else if (ng < MIN_GRIDS && maxlgrids != ng && dir != -1) {
	 dir = 1;
	 refine_alpha /= 1.1;
	 //			printf("Increasing refinement %i %e\n", refine_flag_cnt() + 1, refine_alpha);
	 done = false;
	 } else {
	 done = true;
	 }
	 } while (!done);*/
	rc = use_refine_flags();
	if (rc) {
		enforce_boundaries();
	}
	return rc;
}

void GridNode::inject_from_children() {
	ChildIndex c;
	GridNode* child;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->inject_from_children();
		}
	}
	for (c = 0; c < OCT_NCHILD; c++) {
		child = static_cast<GridNode*> (get_child(c));
		if (child != NULL) {
			const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
			State cstate;
			int k, j, k0, j0, i, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(cstate,i,j,k,i0,j0,k0)
			for (int index = 0; index <= indexer.max_index(); index++) {
				k = indexer.y(index);
				j = indexer.x(index);
				k0 = (bw + k) / 2 + c.get_z() * (GNX / 2 - bw);
				j0 = (bw + j) / 2 + c.get_y() * (GNX / 2 - bw);
				for (i = bw; i < GNX - bw; i += 2) {
					i0 = (bw + i) / 2 + c.get_x() * (GNX / 2 - bw);
					cstate = State(0.0);
					cstate += (*child)(i + 0, j + 0, k + 0);
					cstate += (*child)(i + 1, j + 0, k + 0);
					cstate += (*child)(i + 0, j + 1, k + 0);
					cstate += (*child)(i + 1, j + 1, k + 0);
					cstate += (*child)(i + 0, j + 0, k + 1);
					cstate += (*child)(i + 1, j + 0, k + 1);
					cstate += (*child)(i + 0, j + 1, k + 1);
					cstate += (*child)(i + 1, j + 1, k + 1);
					cstate *= 0.125;
					(*this)(i0, j0, k0) = cstate;
				}
			}
		}
	}
}

void GridNode::adjust_x_flux() {
	ChildIndex ci;
	const GridNode* child_l;
	const GridNode* child_r;
	const GridNode* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->adjust_x_flux();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(0, cj, ck);
				if (l == 0) {
					child_r = static_cast<const GridNode*> (get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const GridNode*> (get_sibling(XL)->get_child(ci));
					i0 = bw;
				} else if (l == 1) {
					child_l = static_cast<const GridNode*> (get_child(ci));
					ci.set_x(1);
					child_r = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const GridNode*> (get_sibling(XU)->get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX - bw;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = GNX - bw;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = bw;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
					Vector<Real, STATE_NF> v;
					int k, j, k0, j0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,k,j,k0,j0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (bw + k) / 2 + ck * (GNX / 2 - bw);
						j0 = (bw + j) / 2 + cj * (GNX / 2 - bw);
						v = +child->get_flux(i, j + 0, k + 0);
						v += child->get_flux(i, j + 1, k + 0);
						v += child->get_flux(i, j + 0, k + 1);
						v += child->get_flux(i, j + 1, k + 1);
						v *= 0.25;
						set_flux(v, i0, j0, k0);
					}
				}
			}
		}
	}

}

void GridNode::adjust_y_flux() {
	ChildIndex ci;
	const GridNode* child_l;
	const GridNode* child_r;
	const GridNode* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->adjust_y_flux();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, 0, ck);
				if (l == 0) {
					child_r = static_cast<const GridNode*> (get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const GridNode*> (get_sibling(YL)->get_child(ci));
					i0 = bw;
				} else if (l == 1) {
					child_l = static_cast<const GridNode*> (get_child(ci));
					ci.set_y(1);
					child_r = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const GridNode*> (get_sibling(YU)->get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX - bw;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = GNX - bw;

				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = bw;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
					Vector<Real, STATE_NF> v;
					int k, j, k0, j0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,k,j,k0,j0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (bw + k) / 2 + ck * (GNX / 2 - bw);
						j0 = (bw + j) / 2 + cj * (GNX / 2 - bw);
						v = +child->get_flux(j + 0, i, k + 0);
						v += child->get_flux(j + 1, i, k + 0);
						v += child->get_flux(j + 0, i, k + 1);
						v += child->get_flux(j + 1, i, k + 1);
						v *= 0.25;
						set_flux(v, j0, i0, k0);
					}
				}
			}
		}
	}

}

void GridNode::adjust_z_flux() {
	ChildIndex ci;
	const GridNode* child_l;
	const GridNode* child_r;
	const GridNode* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<GridNode*> (get_child(i))->adjust_z_flux();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, ck, 0);
				if (l == 0) {
					child_r = static_cast<const GridNode*> (get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const GridNode*> (get_sibling(ZL)->get_child(ci));
					i0 = bw;
				} else if (l == 1) {
					child_l = static_cast<const GridNode*> (get_child(ci));
					ci.set_z(1);
					child_r = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const GridNode*> (get_sibling(ZU)->get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const GridNode*> (get_child(ci));
					i0 = GNX - bw;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = GNX - bw;

				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = bw;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
					Vector<Real, STATE_NF> v;
					int k, j, k0, j0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,k,j,k0,j0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (bw + k) / 2 + ck * (GNX / 2 - bw);
						j0 = (bw + j) / 2 + cj * (GNX / 2 - bw);
						v = +child->get_flux(j + 0, k + 0, i);
						v += child->get_flux(j + 1, k + 0, i);
						v += child->get_flux(j + 0, k + 1, i);
						v += child->get_flux(j + 1, k + 1, i);
						v *= 0.25;
						set_flux(v, j0, k0, i0);
					}
				}
			}
		}
	}

}

int GridNode::get_zone_cnt() const {
	const GridNode* child;
	int zone_cnt = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		child = static_cast<const GridNode*> (get_child(i));
		if (child != NULL) {
			zone_cnt += child->get_zone_cnt();
		} else {
			zone_cnt += GRID_NZONES / OCT_NCHILD;
		}
	}
	return zone_cnt;
}

int GridNode::nvar_output() const {
	return STATE_NF;
}

void GridNode::load_output(grid_output_t* go, int i, int j, int k) const {
	for (int l = 0; l < STATE_NF; l++) {
		go->ele[l][go->ei] = (*this)(i, j, k)[l];
	}
}

const char* GridNode::output_field_names(int i) const {
	return State::field_name(i);
}

void GridNode::output(grid_output_t* ptr) const {
	for (int k = 0; k < GNX + 1; k++) {
		for (int j = 0; j < GNX + 1; j++) {
			for (int i = 0; i < GNX + 1; i++) {
				*(ptr->x) = xf(i);
				*(ptr->y) = yf(j);
				*(ptr->z) = zf(k);
				ptr->x++;
				ptr->y++;
				ptr->z++;
			}
		}
	}
	for (int k = bw; k < GNX - bw; k++) {
		for (int j = bw; j < GNX - bw; j++) {
			for (int i = bw; i < GNX - bw; i++) {
				if (zone_is_refined(i, j, k)) {
					continue;
				}
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 1)) + ptr->pi;
				this->load_output(ptr, i, j, k);
				ptr->ei++;
			}
		}
	}
	ptr->pi += GRID_NNODES;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<const GridNode*> (get_child(i))->output(ptr);
		}
	}
}

bool GridNode::zone_is_refined(int i, int j, int k) const {
	ChildIndex c;
	i = (2 * i) / GNX;
	j = (2 * j) / GNX;
	k = (2 * k) / GNX;
	c.set_index(i, j, k);
	return (get_child(c) != NULL);
}

GridNode::GridNode() :
	Grid(), OctNode() {
	age = -1;
	int k, j, i;
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	for (i = 0; i < OCT_NCHILD; i++) {
		refine[i] = false;
	}
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			E0(i, j, k) = Vector<Real, STATE_NF> (0.0);
		}
	}

}

void GridNode::init() {
	OctNode* sib;
	assert( get_level() == 0 );
	for (int i = 0; i < OCT_NSIB; i++) {
		sib = this->new_phys_bound();
		this->set_sibling(OctFace(i), sib);
		sib->set_sibling(invert(OctFace(i)), this);
	}
}

int GridNode::ngrids() const {
	int sum = 1;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			sum += static_cast<const GridNode*> (get_child(i))->ngrids();
		} else if (refine[i]) {
			sum++;
		}
	}
	return sum;
}

int GridNode::ngrids(int l) const {
	int sum = 0;
	if (l == get_level()) {
		sum += 1;
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				sum += static_cast<const GridNode*> (get_child(i))->ngrids(l);
			}
		}
	}
	return sum;
}

GridNode::~GridNode() {
}

void GridNode::output(const char* prefix, int counter) const {
	assert( get_level() == 0 );
	char filename[32];
	int cnt, nzones, nnodes;
	grid_output_t go;
	float* coords[3];
	DBfile* db;
	DBoptlist* olist;
	char* coordnames[3];
	int shapesize[1];
	int shapecnt[1];
	int shapetype[1];
	int nshapes = 1;
	float ftime = float(get_time() / t_norm);
	//double dtime = double(get_time());
	const int nf = this->nvar_output();
	sprintf(filename, "%s.%i.silo", prefix, counter);
	olist = DBMakeOptlist(1);
	DBAddOption(olist, DBOPT_TIME, &ftime);
	//	DBAddOption(olist, DBOPT_DTIME, &dtime);
	db = DBCreate(filename, DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
	cnt = get_node_cnt();
	nnodes = cnt * GRID_NNODES;
	nzones = get_zone_cnt();
	for (int i = 0; i < 3; i++) {
		coordnames[i] = new char[2];
		coords[i] = reinterpret_cast<float*> (new OReal[nnodes]);
	}
	go.nodelist = new int[8 * nzones];
	go.x = reinterpret_cast<OReal*> (coords[0]);
	go.y = reinterpret_cast<OReal*> (coords[1]);
	go.z = reinterpret_cast<OReal*> (coords[2]);
	go.ele = new OReal*[nf];
	for (int i = 0; i < nf; i++) {
		go.ele[i] = new OReal[nzones];
	}
	go.pi = 0;
	go.ni = 0;
	go.ei = 0;
	GridNode::output(&go);
	shapesize[0] = 8;
	shapecnt[0] = nzones;
	shapetype[0] = DB_ZONETYPE_HEX;
	strcpy(coordnames[0], "x");
	strcpy(coordnames[1], "y");
	strcpy(coordnames[2], "z");
	DBPutZonelist2(db, "zones", nzones, 3, go.nodelist, 8 * nzones, 0, 0, 0, shapetype, shapesize, shapecnt, nshapes,
			olist);
	DBPutUcdmesh(db, "mesh", 3, coordnames, coords, nnodes, nzones, "zones", NULL, DB_OREAL, olist);
	for (int i = 0; i < nf; i++) {
		DBPutUcdvar1(db, output_field_names(i), "mesh", (float*) go.ele[i], nzones, 0, 0, DB_OREAL, DB_ZONECENT, olist);
	}
	DBClose(db);
	char str[80];
	sprintf(str, "gzip -f %s\n", filename);
	system(str);
	delete[] go.nodelist;
	for (int i = 0; i < nf; i++) {
		delete[] go.ele[i];
	}
	delete[] go.ele;
	for (int i = 0; i < 3; i++) {
		delete coords[i];
		delete coordnames[i];
	}
}

Vector<State, 4> GridNode::state_sum() const {
	Vector<State, 4> s0;
	Real h3 = pow(get_dx(), 3);
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	Real a0, a1, a2, a3;
	for (int l = 0; l < STATE_NF; l++) {
		a0 = a1 = a2 = a3 = 0.0;
		Real ds;
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) reduction(+:a1,a2,a3,a0) private(ds,k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = bw; i < GNX - bw; i++) {
				if (!zone_is_refined(i, j, k)) {
					ds = (*this)(i, j, k)[l] * h3;
					a0 += ds;
					a1 += ds * xc(i);
					a2 += ds * yc(j);
					a3 += ds * zc(k);
				}
			}
		}
		s0[0][l] = a0;
		s0[1][l] = a1;
		s0[2][l] = a2;
		s0[3][l] = a3;
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			s0 += static_cast<const GridNode*> (get_child(i))->state_sum();
		}
	}
	/*	if (get_level() == 0) {
	 printf("%e %e %e\n", s0[0][0], s0[1][0] / s0[0][0], s0[2][0] / s0[0][0]);
	 printf("%e %e %e\n", s0[0][1], s0[1][1] / s0[0][1], s0[2][1] / s0[0][1]);
	 }*/
	return s0;
}

State GridNode::state_max() const {
	Vector<Real, STATE_NF> s0;
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	Real smax;
	for (int l = 0; l < STATE_NF; l++) {
		smax = -1.0E+99;
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = bw; i < GNX - bw; i++) {
				if (!zone_is_refined(i, j, k)) {
#pragma omp critical
					smax = max((*this)(i, j, k)[l], smax);
				}
			}
		}
		s0[l] = smax;
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			State s1 = static_cast<const GridNode*> (get_child(i))->state_max();
			for (int l = 0; l < STATE_NF; l++) {
				s0[l] = max(s1[l], s0[l]);
			}
		}
	}
	return s0;
}

State GridNode::state_min() const {
	Vector<Real, STATE_NF> s0;
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	Real smin;
	for (int l = 0; l < STATE_NF; l++) {
		smin = +1.0E+99;
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = bw; i < GNX - bw; i++) {
				if (!zone_is_refined(i, j, k)) {
#pragma omp critical
					smin = min((*this)(i, j, k)[l], smin);
				}
			}
		}
		s0[l] = smin;
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			State s1 = static_cast<const GridNode*> (get_child(i))->state_min();
			for (int l = 0; l < STATE_NF; l++) {
				s0[l] = min(s1[l], s0[l]);
			}
		}
	}
	return s0;
}

