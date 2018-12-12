#include "ConstNumberPerProc3D.hpp"
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#include "mpi_commands.hpp"
#include <mpi.h>
#include "HilbertProcPositions.hpp"
#include "../3D/r3d/Intersection3D.hpp"
#endif

#ifdef RICH_MPI

namespace
{
	Vector3D GetProcCM(Tessellation3D const& tess)
	{
		size_t Ncor = tess.GetPointNo();
		Vector3D res;
		for (size_t i = 0; i < Ncor; ++i)
			res += tess.GetMeshPoint(i);
		res *= 1.0 / static_cast<double>(std::max(static_cast<size_t>(1), Ncor));
		return res;
	}

	std::array<double, 4> DotQ(std::array<double, 4> const& q1, std::array<double, 4> const& q2)
	{
		std::array<double, 4> res{};
		res[3] = (-q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3]);
		res[0] = (q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1]);
		res[1] = (q1[3] * q2[1] - q1[0] * q2[2] + q1[1] * q2[3] + q1[2] * q2[0]);
		res[2] = (q1[3] * q2[2] + q1[0] * q2[1] - q1[1] * q2[0] + q1[2] * q2[3]);
		return res;
	}

	std::array<double, 9> GetMatrix(std::array<double, 4> const& q)
	{
		std::array<double, 9> res = { 1 - 2 * q[1] * q[1] - 2 * q[2] * q[2], 2 * q[0] * q[1] - 2 * q[2] * q[3], 2 * q[0] * q[2] + 2 * q[1] * q[3],
			2 * q[0] * q[1] + 2 * q[2] * q[3], 1 - 2 * q[0] * q[0] - 2 * q[2] * q[2], 2 * q[1] * q[2] - 2 * q[0] * q[3],
			2 * q[0] * q[2] - 2 * q[1] * q[3], 2 * q[1] * q[2] + 2 * q[0] * q[3], 1 - 2 * q[0] * q[0] - 2 * q[1] * q[1] };
		return res;
	}

	std::array<double, 9> MatrixMult(std::array<double, 9> const&Q, std::array<double, 9> const& A)
	{
		std::array<double, 9> res{};
		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				for (size_t k = 0; k < 3; ++k)
					res[i * 3 + j] += Q[i * 3 + k] * A[j + k * 3];
			}
		}
		return res;
	}

	std::array<double, 9> MatrixMultT(std::array<double, 9> const&Q, std::array<double, 9> const& A)
	{
		std::array<double, 9> res{};
		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				res[i * 3 + j] = A[j * 3 + i];
			}
		}
		return MatrixMult(Q, res);
	}

	void eigsrt(std::array<double, 3> &d, std::array<double, 9> &v)
	{
		int k;
		int n = d.size();
		for (int i = 0; i < n - 1; i++) 
		{
			double p = d[k = i];
			for (int j = i; j < n; j++)
				if (d[j] >= p) p = d[k = j];
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (int j = 0; j < n; j++) 
				{
					p = v[j*3 + i];
					v[j*3 + i] = v[j*3 + k];
					v[j*3 + k] = p;
				}
			}
		}
	}

	struct Jacobi 
	{
		const int n;
		std::array<double, 9> a, v;
		std::array<double, 3> d;
		int nrot;
		const double EPS;

		Jacobi(std::array<double, 9> const& aa) : n(3), a(aa), v(), d(), nrot(0),
			EPS(std::numeric_limits<double>::epsilon())
		{
			int i, j, ip, iq;
			double tresh, theta, tau, t, sm, s, h, g, c;
			std::array<double, 3> b,z;
			for (ip = 0; ip < n; ip++) 
			{
				for (iq = 0; iq < n; iq++) 
					v[ip*3+ iq] = 0.0;
				v[ip*3+ip] = 1.0;
			}
			for (ip = 0; ip < n; ip++) 
			{
				b[ip] = d[ip] = a[ip*3 + ip];
				z[ip] = 0.0;
			}
			for (i = 1; i <= 50; i++) 
			{
				sm = 0.0;
				for (ip = 0; ip < n - 1; ip++) {
					for (iq = ip + 1; iq < n; iq++)
						sm += abs(a[ip*3 + iq]);
				}
				if (sm == 0.0) 
				{
					eigsrt(d, v);
					return;
				}
				if (i < 4)
					tresh = 0.2*sm / (n*n);
				else
					tresh = 0.0;
				for (ip = 0; ip < n - 1; ip++) 
				{
					for (iq = ip + 1; iq < n; iq++) 
					{
						g = 100.0*abs(a[ip*3 + iq]);
						if (i > 4 && g <= EPS * abs(d[ip]) && g <= EPS * abs(d[iq]))
							a[ip*3 + iq] = 0.0;
						else 
							if (abs(a[ip*3 + iq]) > tresh) 
						{
							h = d[iq] - d[ip];
							if (g <= EPS * abs(h))
								t = (a[ip*3 + iq]) / h;
							else 
							{
								theta = 0.5*h / (a[ip*3 + iq]);
								t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
								if (theta < 0.0) t = -t;
							}
							c = 1.0 / sqrt(1 + t * t);
							s = t * c;
							tau = s / (1.0 + c);
							h = t * a[ip*3 + iq];
							z[ip] -= h;
							z[iq] += h;
							d[ip] -= h;
							d[iq] += h;
							a[ip*3 + iq] = 0.0;
							for (j = 0; j < ip; j++)
								rot(a, s, tau, j, ip, j, iq);
							for (j = ip + 1; j < iq; j++)
								rot(a, s, tau, ip, j, j, iq);
							for (j = iq + 1; j < n; j++)
								rot(a, s, tau, ip, j, iq, j);
							for (j = 0; j < n; j++)
								rot(v, s, tau, j, ip, j, iq);
							++nrot;
						}
					}
				}
				for (ip = 0; ip < n; ip++) 
				{
					b[ip] += z[ip];
					d[ip] = b[ip];
					z[ip] = 0.0;
				}
			}
			throw("Too many iterations in routine jacobi");
		}

		inline void rot(std::array<double, 9> &a, const double s, const double tau, const int i,
			const int j, const int k, const int l)
		{
			double g = a[i*3 + j];
			double h = a[k*3 + l];
			a[i*3 + j] = g - s * (h + g * tau);
			a[k*3 + l] = h + s * (g - h * tau);
		}
	};

	std::array<double, 4> Diagonalizer(std::array<double, 9> const &A, std::array<double, 3> &eigv)
	{
		// A must be a symmetric matrix.
		// returns quaternion q such that its corresponding matrix Q 
		// can be used to Diagonalize A
		// Diagonal matrix D = Q * A * Transpose(Q);  and  A = QT*D*Q
		// The rows of q are the eigenvectors D's diagonal is the eigenvalues
		// As per 'row' convention if double3x3 Q = q.getmatrix(); then v*Q = q*v*conj(q)
		int maxsteps = 240;  // certainly wont need that many.
		int i;
		std::array<double, 4> q{ 0, 0, 0, 1 };
		for (i = 0; i < maxsteps; i++)
		{
			std::array<double, 9> Q = GetMatrix(q); // v*Q == q*v*conj(q)
			std::array<double, 9> D = MatrixMult(Q, MatrixMultT(A, Q));
			eigv[0] = D[0];
			eigv[1] = D[4];
			eigv[2] = D[8];
			std::array<double, 3> offdiag = { D[5], D[2], D[1] }; // elements not on the diagonal
			std::array<double, 3> om = { std::fabs(offdiag[0]), std::fabs(offdiag[1]), std::fabs(offdiag[2]) }; // mag of each offdiag elem
			int k = (om[0] > om[1] && om[0] > om[2]) ? 0 : (om[1] > om[2]) ? 1 : 2; // index of largest element of offdiag
			int k1 = (k + 1) % 3;
			int k2 = (k + 2) % 3;
			if (offdiag[k] == 0.0)
				break;  // diagonal already
			double thet = (D[k2 * 3 + k2] - D[k1 * 3 + k1]) / (2.0*offdiag[k]);
			double sgn = (thet > 0.0) ? 1.0 : -1.0;
			thet *= sgn; // make it positive
			double t = sgn / (thet + ((thet < 1e6) ? std::sqrt(thet*thet + 1.0) : thet + 0.5 / thet)); // sign(T)/(|T|+sqrt(T^2+1))
			double c = 1.0 / std::sqrt(t*t + 1.0); //  c= 1/(t^2+1) , t=s/c 
			if (c > (1 - 1e-13))
				break;  // no room for improvement - reached machine precision.
			std::array<double, 4> jr = { 0, 0, 0, 0 }; // jacobi rotation for this iteration.
			if (std::abs(t) < 1e-6)
				jr[k] = sgn * std::abs(t*0.5);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
			else
				jr[k] = sgn * std::sqrt((1.0 - c) / 2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)  
			jr[3] = std::abs(jr[k]) < 1e-6 ? 1 - 0.5*jr[k] * jr[k] : std::sqrt(1.0 - jr[k] * jr[k]);
			if (jr[3] > (1 - 1e-13))
				break; // reached limits of doubleing point precision
			q = DotQ(q, jr);
			double sum = 0;
			for (size_t j = 0; j < 4; ++j)
				sum += q[j] * q[j];
			sum = 1.0 / std::sqrt(sum);
			for (size_t j = 0; j < 4; ++j)
				q[j] *= sum;
		}
		assert(i < maxsteps);
		return q;
	}

	std::array<double, 9> GetMomentOfInertia(Tessellation3D const& tess, Vector3D const& CM)
	{
		std::array<double, 9> res{};
		size_t N = tess.GetPointNo();
		for (size_t i = 0; i < N; ++i)
		{
			Vector3D p = tess.GetMeshPoint(i) - CM;
			res[0] += p.y*p.y + p.z*p.z;
			res[1] -= p.x*p.y;
			res[2] -= p.x*p.z;
			res[4] += p.x*p.x + p.z*p.z;
			res[5] -= p.z*p.y;
			res[8] += p.y*p.y + p.x*p.x;
		}
		res[3] = res[1];
		res[6] = res[2];
		res[7] = res[5];
		return res;
	}

	Vector3D GetMinEigVector(Tessellation3D const& tess, Vector3D const& CM)
	{
		std::array<double, 9> I = GetMomentOfInertia(tess, CM);
		std::array<double, 3> eigv;
		//std::array<double, 4> q = Diagonalizer(I, eigv);
		//std::array<double, 9> Q = GetMatrix(q);
		Jacobi jacobi(I);
		std::array<double, 9> Q2 = jacobi.v;
	//	size_t index = std::min_element(eigv.begin(), eigv.end()) - eigv.begin();
		//Vector3D res(Q[index * 3], Q[index * 3 + 1], Q[index * 3 + 2]);
		Vector3D res(Q2[2], Q2[5], Q2[8]);
		return res;
	}

	Vector3D CalcGradient(Tessellation3D const& tess, Tessellation3D const& tproc, vector<int> const& NPerProc)
	{
		Vector3D res;
		size_t N = tess.GetPointNo();
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		double totalV = tproc.GetVolume(static_cast<size_t>(rank));
		// Do default result by pointing atneihbor with largest N
		int MaxN = 0;
		size_t loc = 0;
		std::vector<size_t> neigh = tproc.GetNeighbors(static_cast<size_t>(rank));
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (neigh[i] < tproc.GetPointNo())
			{
				if (NPerProc[neigh[i]] > MaxN)
				{
					MaxN = NPerProc[neigh[i]];
					loc = neigh[i];
				}
			}
		}
		res = (MaxN / tproc.GetVolume(loc) - N / totalV)*(tproc.GetCellCM(loc) - tproc.GetCellCM(static_cast<size_t>(rank))) /
			ScalarProd(tproc.GetCellCM(loc) - tproc.GetCellCM(static_cast<size_t>(rank)), tproc.GetCellCM(loc) - tproc.GetCellCM(static_cast<size_t>(rank)));
		// deal with 0 points
		if (N == 0)
			return res;

		Vector3D CM = GetProcCM(tess);
		Vector3D eigv = GetMinEigVector(tess, CM);
		Vector3D up, down;

		size_t counter = 0;
		for (size_t i = 0; i < N; ++i)
		{
			Vector3D p = tess.GetMeshPoint(i) - CM;
			if (ScalarProd(p, eigv) > 0)
				down += p;
			else
			{
				up += p;
				++counter;
			}
		}
		if (counter == 0 || counter == N)
			return res;

		// Find volumes of up and down
		r3d_plane plane;
		plane.n.xyz[0] = eigv.x;
		plane.n.xyz[1] = eigv.y;
		plane.n.xyz[2] = eigv.z;
		plane.d = -1 * ScalarProd(eigv, CM);
		r3d_poly poly;
		std::vector<size_t> itemp, all_indeces;
		std::vector<std::vector<int> > faceinds;

		if (GetPoly(tproc, static_cast<size_t>(rank), poly, itemp, all_indeces, faceinds))
		{
			r3d_clip(&poly, &plane, 1);
			if (poly.nverts > 0)
			{
				double m[1];
				m[0] = 0;
				r3d_reduce(&poly, m, 0);
				double Vup = std::abs(m[0]);
				double dn = counter / Vup - (N - counter) / (totalV - Vup);
				Vector3D dr = up / counter - down / (N - counter);
				res = dn * dr / ScalarProd(dr, dr);
			}
			return res;
		}
		return res;
	}
}
#endif

ConstNumberPerProc3D::~ConstNumberPerProc3D(void) {}


ConstNumberPerProc3D::ConstNumberPerProc3D(double speed, double RoundSpeed, int mode, bool Hilbert) :
	speed_(speed), RoundSpeed_(RoundSpeed), mode_(mode), Hilbert_(Hilbert), run_counter_(100) {}

#ifdef RICH_MPI
void ConstNumberPerProc3D::Update(Tessellation3D& tproc, Tessellation3D const& tlocal)const
{
	int nproc = static_cast<int>(tproc.GetPointNo());
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<double> R(static_cast<size_t>(nproc));
	double dx = 0;
	double dy = 0;
	double dz = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		R[i] = tproc.GetWidth(i);
	// Make cell rounder
	const Vector3D& CM = tproc.GetCellCM(rank);
	Vector3D point = tproc.GetMeshPoint(rank);
	// Find out how many points each proc has
	vector<int> NPerProc(static_cast<size_t>(nproc));
	int mypointnumber = static_cast<int>(tlocal.GetPointNo() + 1);

	MPI_Gather(&mypointnumber, 1, MPI_INT, &NPerProc[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&NPerProc[0], nproc, MPI_INT, 0, MPI_COMM_WORLD);

	int IdealPerProc = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		IdealPerProc += NPerProc[i];
	IdealPerProc /= nproc;
	double load = 0;
	for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		load = std::max(load, static_cast<double>(NPerProc[i]) / static_cast<double>(IdealPerProc));

	double NewSpeed = speed_;
	if (load > 3)
		NewSpeed *= std::min(std::pow(load - 2, 3), 0.15 / speed_);

	const double d = fastabs(CM - tproc.GetMeshPoint(rank));
	double dxround = 0, dyround = 0, dzround = 0;
	if (d > 0.1*R[static_cast<size_t>(rank)])
	{
		dxround = RoundSpeed_ * NewSpeed*(CM.x - point.x);
		dyround = RoundSpeed_ * NewSpeed*(CM.y - point.y);
		dzround = RoundSpeed_ * NewSpeed*(CM.z - point.z);
	}

	if (Hilbert_ && load > 3.75 && (run_counter_ % 100 == 0))
	{
		vector<Vector3D> res = HilbertProcPositions(tlocal);
		tproc.Build(res);
		return;
	}
	++run_counter_;
	Vector3D RankCM = GetProcCM(tlocal);
	//Vector3D RankCM = tproc.GetCellCM(static_cast<size_t>(rank));
	vector<double> tosend = list_serialize(vector<Vector3D>(1, RankCM));
	vector<double> torecv(static_cast<size_t>(nproc) * 3);
	MPI_Gather(&tosend[0], 3, MPI_DOUBLE, &torecv[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector3D> RankCMs = list_unserialize(torecv, RankCM);
	double MyR = R[static_cast<size_t>(rank)];
	// Move point according to density
	if (mode_ == 1 || mode_ == 3)
	{
		for (size_t i = 0; i < static_cast<size_t>(nproc); ++i)
		{
			if (i == static_cast<size_t>(rank))
				continue;
			Vector3D otherpoint = RankCMs[i];
			double dist = std::sqrt((point.x - otherpoint.x)*(point.x - otherpoint.x) +
				(point.y - otherpoint.y)*(point.y - otherpoint.y) + (point.z - otherpoint.z)*(point.z - otherpoint.z)
				+ 0.5*MyR * R[i]);
			double min_volume = 0.35*std::min(tproc.GetVolume(i), 0.5*dist*dist*dist);
			double merit = static_cast<double>(NPerProc[i] - IdealPerProc) / IdealPerProc;
			merit = std::max(merit, -0.5);
			double temp = merit * (point.x - otherpoint.x) * min_volume / (dist*dist*dist);
			dx -= NewSpeed * temp;
			temp = merit * (point.y - otherpoint.y) * min_volume / (dist*dist*dist);
			dy -= NewSpeed * temp;
			temp = merit * (point.z - otherpoint.z) * min_volume / (dist*dist*dist);
			dz -= NewSpeed * temp;
		}
	}
	point = tproc.GetMeshPoint(static_cast<size_t>(rank));
	double old_dx = dx;
	double old_dy = dy;
	double old_dz = dz;
	dx = 0;
	dy = 0;
	dz = 0;
	// Moving according to pressure
	if (mode_ == 1 || mode_ == 2)
	{
		vector<size_t> neigh = tproc.GetNeighbors(rank);
		size_t Nneigh = neigh.size();

		const double neigheps = 0.05;
		for (size_t i = 0; i < Nneigh; ++i)
		{
			if (static_cast<int>(neigh[i]) >= nproc)
				continue;
			const double dist = fastabs(point - tproc.GetMeshPoint(neigh[i]));
			const double mind = neigheps * std::min(MyR, R[neigh[i]]);
			if (dist < mind)
			{
				double speed2 = NewSpeed * 10 * (mind - dist) *MyR / (dist*dist);
				dx += speed2 * (point.x - tproc.GetMeshPoint(neigh[i]).x);
				dy += speed2 * (point.y - tproc.GetMeshPoint(neigh[i]).y);
				dz += speed2 * (point.z - tproc.GetMeshPoint(neigh[i]).z);
			}
			else
			{
				Vector3D otherpoint = RankCMs[neigh[i]];
				double merit = static_cast<double>(NPerProc[static_cast<size_t>(rank)] - NPerProc[neigh[i]]) / IdealPerProc;
				double dr = fastabs(otherpoint - point);
				dx -= NewSpeed * merit*(otherpoint.x - point.x)*MyR / dr;
				dy -= NewSpeed * merit*(otherpoint.y - point.y)*MyR / dr;
				dz -= NewSpeed * merit*(otherpoint.z - point.z)*MyR / dr;
			}
		}
	}

	if (mode_ == 4)
	{
		dx = 0;
		dy = 0;
		dz = 0;
		Vector3D grad = CalcGradient(tlocal, tproc, NPerProc);
		double prefactor = (mypointnumber - IdealPerProc) / (ScalarProd(grad, grad)*MyR*MyR + 1e-8*IdealPerProc / MyR);
		dx = prefactor * grad.x;
		dy = prefactor * grad.y;
		dz = prefactor * grad.z;
	}
	old_dx += dx;
	old_dy += dy;
	old_dz += dz;

	// Get min distance to other points
	vector<size_t> neigh = tproc.GetNeighbors(rank);
	size_t Nneigh = neigh.size();
	double mind_1 = 0;
	for (size_t i = 0; i < Nneigh; ++i)
	{
		if (static_cast<int>(neigh[i]) >= nproc)
			continue;
		const double dist = fastabs(point - tproc.GetMeshPoint(neigh[i]));
		mind_1 = std::max(mind_1, 1.0 / dist);
	}
	double maxR = std::min(MyR, 1.5 / mind_1);
	double r_dx = fastabs(Vector3D(old_dx, old_dy, old_dz));
	if (r_dx > NewSpeed*maxR)
	{
		old_dx *= NewSpeed * maxR / r_dx;
		old_dy *= NewSpeed * maxR / r_dx;
		old_dz *= NewSpeed * maxR / r_dx;
	}
	// Add the round cells
	double round_reduce = std::min(1.0, maxR / MyR);
	old_dx += dxround * round_reduce;
	old_dy += dyround * round_reduce;
	old_dz += dzround * round_reduce;
	r_dx = fastabs(Vector3D(old_dx, old_dy, old_dz));
	if (r_dx > NewSpeed*maxR)
	{
		old_dx *= NewSpeed * maxR / r_dx;
		old_dy *= NewSpeed * maxR / r_dx;
		old_dz *= NewSpeed * maxR / r_dx;
	}
	dx = old_dx;
	dy = old_dy;
	dz = old_dz;
	// Make sure not out of bounds
	point = tproc.GetMeshPoint(rank);
	const double close = 0.999;
	const double wx = tproc.GetWidth(rank);
	Vector3D ll = tproc.GetBoxCoordinates().first;
	Vector3D ur = tproc.GetBoxCoordinates().second;
	const Vector3D center(0.5*(ll + ur));
	if (point.x + dx > (ur.x - (1 - close)*wx))
	{
		if ((ur.x - point.x) < ((1 - close)*wx))
			dx = -wx * (1 - close);
		else
			dx = 0.5*(ur.x - point.x);
	}
	if ((point.x + dx) < (ll.x + (1 - close)*wx))
	{
		if ((-ll.x + point.x) < ((1 - close)*wx))
			dx = wx * (1 - close);
		else
			dx = -0.5*(-ll.x + point.x);
	}
	if (point.y + dy > (ur.y - (1 - close)*wx))
	{
		if ((ur.y - point.y) < ((1 - close)*wx))
			dy = -wx * (1 - close);
		else
			dy = 0.5*(ur.y - point.y);
	}
	if ((point.y + dy) < (ll.y + (1 - close)*wx))
	{
		if ((-ll.y + point.y) < ((1 - close)*wx))
			dy = wx * (1 - close);
		else
			dy = -0.5*(-ll.y + point.y);
	}
	if (point.z + dz > (ur.z - (1 - close)*wx))
	{
		if ((ur.z - point.z) < ((1 - close)*wx))
			dz = -wx * (1 - close);
		else
			dz = 0.5*(ur.z - point.z);
	}
	if ((point.z + dz) < (ll.z + (1 - close)*wx))
	{
		if ((-ll.z + point.z) < ((1 - close)*wx))
			dz = wx * (1 - close);
		else
			dz = -0.5*(-ll.z + point.z);
	}
	Vector3D cor = tproc.GetMeshPoint(static_cast<size_t>(rank)) + Vector3D(dx, dy, dz);
	// Have all processors have the same points
	tosend = list_serialize(vector<Vector3D>(1, cor));
	MPI_Gather(&tosend[0], 3, MPI_DOUBLE, &torecv[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&torecv[0], nproc * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<Vector3D> cortemp = list_unserialize(torecv, cor);
	tproc.Build(cortemp);
}
#endif
