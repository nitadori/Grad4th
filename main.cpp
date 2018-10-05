#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include "Nbody.h"

int main(int argc, char **argv){
	Nbody sys;

	assert(argc > 1);
	FILE *fp = fopen(argv[1], "r");
	assert(fp);
	sys.read(fp);
	fclose(fp);

	double en0 = sys.energy(stderr);

	sys.calc_acc();

	const double tick = 1.0 / 8192.0;

	double err_max = 0.0;

	for(int i=0; i<100; i++){
		// sys.DKD_2nd(tick);
		// sys.KDK_2nd(tick);
		// sys.KDKDK_2nd(tick);
		sys.KDKDK_4th(tick);
		double en1 = sys.energy();
		double de = (en1 - en0) / en0;
		printf("%e %e\n", sys.tsys, de);

		err_max = std::max(err_max, fabs(de));
	}
	printf("error max: %e\n", err_max);

	return 0;
}


