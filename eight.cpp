#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#define SKIP_SELF
#include "Nbody.h"

#ifndef INTGRT
#define INTGRT KDKDK_4th
#endif

int main(int argc, char **argv){
	Nbody sys;
	sys.resize(3);

	assert(argc > 1);
	const int invtick = atoi(argv[1]);
	const double tick = 1.0 / invtick;

	sys.set_eps(0.0);
	sys.ptcl[0].pos = { 0.97000436, -0.24308753, 0.0};
	sys.ptcl[1].pos = {-0.97000436, +0.24308753, 0.0};
	sys.ptcl[2].pos = {0.0,          0.0,        0.0};
	sys.ptcl[0].vel = {-0.93240737/(-2.0),-0.86473146/(-2.0),0.0};
	sys.ptcl[1].vel = {-0.93240737/(-2.0),-0.86473146/(-2.0),0.0};
	sys.ptcl[2].vel = {-0.93240737,       -0.86473146,       0.0};
	sys.ptcl[0].mass = 1.0;
	sys.ptcl[1].mass = 1.0;
	sys.ptcl[2].mass = 1.0;
	sys.ptcl[0].id = 0;
	sys.ptcl[1].id = 1;
	sys.ptcl[2].id = 2;

	sys.calc_acc();
	double en0 = sys.energy(stderr);
	double err_max = 0.0;

	FILE *fp = fopen("orbit.dat", "w");

	while(sys.tsys < 100.0){
		sys.INTGRT(tick);
		double en1 = sys.energy();
		double de = (en1 - en0) / en0;
		printf("%e %e\n", sys.tsys, de);

#if 1
		fprintf(fp, "%e %e %f %f %f %f %f %f\n",
				sys.tsys, de, 
				sys.ptcl[0].pos.x, sys.ptcl[0].pos.y,
				sys.ptcl[1].pos.x, sys.ptcl[1].pos.y,
				sys.ptcl[2].pos.x, sys.ptcl[2].pos.y);
#endif

		err_max = std::max(err_max, fabs(de));
	}
	printf("%e %e (dt,error_max)\n", 12.*tick, err_max);

	fclose(fp);
}
