#include <cstdio>
#include <cmath>

struct Vec3{
	double x, y, z;

	friend Vec3 operator*(const double s, const Vec3 v){
		return {s*v.x, s*v.y, s*v.z};
	}
	Vec3 &operator+=(const Vec3 rhs){ 
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;

		return *this;
	}

	Vec3 operator-(const Vec3 rhs) const { 
		return {
			x - rhs.x,
			y - rhs.y,
			z - rhs.z
		};
	}

	Vec3 operator+(const Vec3 rhs) const { 
		return {
			x + rhs.x,
			y + rhs.y,
			z + rhs.z
		};
	}

	double operator*(const Vec3 rhs) const {
		return x*rhs.x + y*rhs.y + z*rhs.z;
	}
};

struct Particle{
	long   id;
	double mass;
	double pot;
	Vec3 pos;
	Vec3 vel;
	Vec3 acc;
	Vec3 acorr; // a \dot grad a
	// 15 DP, 120 byte

	void kick(const double h){
		vel += h * acc;
	}

	void drift(const double h){
		pos += h * vel;
	}

	void kick_with_corr_inner(const double h, const double hh){
		vel += h * (acc + (1.0 * hh) * acorr);
		if(id == -1){
			auto ac = hh * acorr;
			printf ("(%e, %e, %e) + (%e, %e, %e)\n",
					acc.x, acc.y, acc.z,
					ac.x, ac.y, ac.z);
		}
	}

	void kick_with_corr(const double tick){
		// tick      =  h / 12;
		// 2h/3      =  8 * tick;
		// h^2 / 48  =  3 * tick * tick;
		kick_with_corr_inner(8.0 * tick, 3.0 * tick * tick);
	}

	void calc_acc(const Particle *pj, const int n, const double eps2){
		Vec3 acc = {0.0, 0.0, 0.0};
		for(int j=0; j<n; j++){
#ifdef SKIP_SELF
			if(id == pj[j].id) continue;
#endif
			Vec3 dr = pj[j].pos - pos;
			double r2 = eps2 + dr*dr;
			double ri2 = 1.0 / r2;
			double ri1 = sqrt(ri2);
			double mri1 = pj[j].mass * ri1;
			double mri3 = mri1 * ri2;

			acc += mri3 * dr;
		}
		this->acc = acc;
	}

	void calc_acorr(const Particle *pj, const int n, const double eps2){
		Vec3 acorr = {0.0, 0.0, 0.0};
		for(int j=0; j<n; j++){
#ifdef SKIP_SELF
			if(id == pj[j].id) continue;
#endif
			Vec3 dr = pj[j].pos - pos;
			Vec3 da = pj[j].acc - acc;
			
			double r2 = eps2 + dr*dr;
			double drda = dr * da;

			double ri2 = 1.0 / r2;
			double ri1 = sqrt(ri2);
			double mri1 = pj[j].mass * ri1;
			double mri3 = mri1 * ri2;

			double alpha = 3.0 * drda * ri2;

			acorr += mri3 * (da - alpha * dr);
		}
		this->acorr = 2.0 * acorr;
	}

	void calc_pot(const Particle *pj, const int n, const double eps2){
		double pot = 0.0;
		for(int j=0; j<n; j++){
			Vec3 dr = pj[j].pos - pos;
			double r2 = eps2 + dr*dr;
			double ri2 = 1.0 / r2;
			double ri1 = sqrt(ri2);
			double mri1 = pj[j].mass * ri1;

			pot -= (id != pj[j].id) ? mri1 : 0.0;
		}
		this->pot = pot;
	}
		 
};
