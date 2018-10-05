#include <cstdio>
#include <cassert>

#include "Particle.h"

struct Nbody{
	int snapid;
	int nbody;
	double tsys;

	Particle *ptcl;

	double eps2;

	void read(FILE *fp){
		assert( 3 == fscanf(fp, "%d %d %lf", &snapid,  &nbody, &tsys) );
		ptcl = new Particle[nbody];

		for(int i=0; i<nbody; i++){
			Particle &p = ptcl[i];
			assert(9 == fscanf(
						fp,
						"%ld %lf %lA %lA %lA %lA %lA %lA %lA",
						&p.id, &p.mass, 
						&p.pos.x, &p.pos.y, &p.pos.z, 
						&p.vel.x, &p.vel.y, &p.vel.z,
						&p.pot) );

		}

		// set default softening here
		double eps = 1.0 / 64.0;
		this->eps2 = eps * eps;

		fprintf(stderr, "read snapshot, N=%d\n", nbody);
	}

	double energy(FILE *fp = NULL){
		for(int i=0; i<nbody; i++){
			ptcl[i].calc_pot(ptcl, nbody, eps2);
		}

		double ke = 0.0, pe = 0.0;
		for(int i=0; i<nbody; i++){
			ke += ptcl[i].mass * (ptcl[i].vel * ptcl[i].vel);
			pe += ptcl[i].mass * ptcl[i].pot;
		}
		ke *= 0.5;
		pe *= 0.5;

		if(fp){
			fprintf(fp, "t = %f, ke = %f, pe = %f, ke+pe = %f\n", tsys, ke, pe, ke+pe);
		}

		return ke + pe;
	}

	void calc_acc(){
#pragma omp parallel for
		for(int i=0; i<nbody; i++){
			ptcl[i].calc_acc(ptcl, nbody, eps2);
		}
	}
	void calc_acorr(){
#pragma omp parallel for
		for(int i=0; i<nbody; i++){
			ptcl[i].calc_acorr(ptcl, nbody, eps2);
		}
	}

	void kick(double h){ 
		for(int i=0; i<nbody; i++) ptcl[i].kick(h);
	}

	void drift(double h){ 
		for(int i=0; i<nbody; i++) ptcl[i].drift(h);
	}

	void DKD_2nd(const double tick){
		double h  = 6.0 * tick;
		double h2 = 3.0 * tick;

		drift(h2);

		calc_acc();
		kick (h );

		drift(h );

		calc_acc();
		kick (h );

		drift(h2);

		tsys += 12.0 * tick;
	}
	void KDK_2nd(const double tick){
		double h  = 6.0 * tick;
		double h2 = 3.0 * tick;

		kick (h2);
		drift(h );

		calc_acc();
		kick (h );

		drift(h );

		calc_acc();
		kick (h2);

		tsys += 12.0 * tick;
	}
	void KDKDK_2nd(const double tick){
		kick (2.0 * tick);
		drift(6.0 * tick);

		calc_acc();
		kick (8.0 * tick);
		drift(6.0 * tick);

		calc_acc();
		kick (2.0 * tick);

		tsys += 12.0 * tick;
	}

	void KDKDK_4th(const double tick){
		kick (2.0 * tick);
		drift(6.0 * tick);

		calc_acc();
		calc_acorr();
		for(int i=0; i<nbody; i++) ptcl[i].kick_with_corr(tick);

		drift(6.0 * tick);

		calc_acc();
		kick (2.0 * tick);

		tsys += 12.0 * tick;
	}
};
