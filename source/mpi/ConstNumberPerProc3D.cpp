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
	/*Vector3D GetNeighborGrad(Tessellation3D const& tproc, std::vector<int> const& Nproc,size_t rank)
	{
		Vector3D grad_temp;
		std::vector<size_t> neigh = tproc.GetNeighbors(rank);
		std::vector<size_t> faces = tproc.GetCellFaces(rank);
		size_t n = neigh.size();
		Vector3D cell_cm = tproc.GetCellCM(rank);
		Vector3D center = tproc.GetMeshPoint(rank);
		std::vector<Vector3D> neighbor_centers(n), neigh_cm(n);
		std::vector<size_t> neighbors(n);
		size_t cell = Nproc[rank];
		for (size_t i = 0; i < n;++i)
		{
			if (neigh[i] < tproc.GetPointNo())
				neighbors[i] = static_cast<size_t>(Nproc[neigh[i]]);
			else
				neighbors[i] = static_cast<size_t>(Nproc[rank]);
			neigh_cm[i] = tproc.GetCellCM(neigh[i]);
			neighbor_centers[i] = tproc.GetMeshPoint(neigh[i]);
		}
		std::vector<Vector3D> c_ij;
		// Create the matrix to invert and the vector to compare
		std::array<double, 9>  m;
		std::fill_n(m.begin(), 9, 0.0);
		c_ij.resize(n);

		for (size_t i = 0; i < n; i++)
		{
			c_ij[i] = neigh_cm[i];
			c_ij[i] += cell_cm;
			c_ij[i] *= -0.5;
			c_ij[i] += tproc.FaceCM(faces[i]);
			const Vector3D r_ij = normalize(neighbor_centers[i] - center);
			const double A = tproc.GetArea(faces[i]);
			m[0] -= c_ij[i].x*r_ij.x*A;
			m[1] -= c_ij[i].y*r_ij.x*A;
			m[2] -= c_ij[i].z*r_ij.x*A;
			m[3] -= c_ij[i].x*r_ij.y*A;
			m[4] -= c_ij[i].y*r_ij.y*A;
			m[5] -= c_ij[i].z*r_ij.y*A;
			m[6] -= c_ij[i].x*r_ij.z*A;
			m[7] -= c_ij[i].y*r_ij.z*A;
			m[8] -= c_ij[i].z*r_ij.z*A;
			if (i == 0)
			{
				grad_temp.x = neighbors[i];
				grad_temp.x *= r_ij.x*A;
				grad_temp.y = neighbors[i];
				grad_temp.y *= r_ij.y*A;
				grad_temp.z = neighbors[i];
				grad_temp.z *= r_ij.z*A;
			}
			else
			{
				grad_temp.x += neighbors[i] * r_ij.x*A;
				grad_temp.y += neighbors[i] * r_ij.y*A;
				grad_temp.z += neighbors[i] * r_ij.z*A;
			}
			grad_temp.x += cell * r_ij.x*A;
			grad_temp.y += cell * r_ij.y*A;
			grad_temp.z += cell * r_ij.z*A;
		}
		double v_inv = 1.0 / tproc.GetVolume(rank);
		for (size_t i = 0; i < 9; ++i)
			m[i] *= v_inv;
		m[0] += 1;
		m[4] += 1;
		m[8] += 1;
		// Find the det
		const double det = -m[2] * m[4] * m[6] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[0] * m[5] * m[7] -
			m[1] * m[3] * m[8] + m[0] * m[4] * m[8];
		// Check none singular
		if (std::abs(det) < 1e-10)
		{
			UniversalError eo("Singular matrix");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			eo.AddEntry("Cell z cor", center.z);
			eo.AddEntry("Cell CMx cor", cell_cm.x);
			eo.AddEntry("Cell CMy cor", cell_cm.y);
			eo.AddEntry("Cell CMz cor", cell_cm.z);
			eo.AddEntry("Cell volume", tproc.GetVolume(rank));
			eo.AddEntry("Det was", det);
			for (size_t i = 0; i < faces.size(); ++i)
			{
				c_ij[0] = tproc.FaceCM(faces[i]) - 0.5 * (neigh_cm[i] + cell_cm);
				eo.AddEntry("Neighbor x", neighbor_centers[i].x);
				eo.AddEntry("Neighbor y", neighbor_centers[i].y);
				eo.AddEntry("Neighbor z", neighbor_centers[i].z);
				eo.AddEntry("Face", static_cast<double>(faces[i]));
				eo.AddEntry("Neighbor Cx", c_ij[0].x);
				eo.AddEntry("Neighbor Cy", c_ij[0].y);
				eo.AddEntry("Neighbor Cz", c_ij[0].z);
				eo.AddEntry("Face Cx", tproc.FaceCM(faces[i]).x);
				eo.AddEntry("Face Cy", tproc.FaceCM(faces[i]).y);
				eo.AddEntry("Face Cz", tproc.FaceCM(faces[i]).z);
				eo.AddEntry("Face area", tproc.GetArea(faces[i]));
			}
			for (size_t i = 0; i < 9; ++i)
				eo.AddEntry("M", m[i]);
			throw eo;
		}
		// Invert the matrix
		std::array<double, 9>  m_inv;
		std::fill_n(m_inv.begin(), 9, 0);
		m_inv[0] = m[4] * m[8] - m[5] * m[7];
		m_inv[1] = m[2] * m[7] - m[1] * m[8];
		m_inv[2] = m[1] * m[5] - m[2] * m[4];
		m_inv[3] = m[5] * m[6] - m[3] * m[8];
		m_inv[4] = m[0] * m[8] - m[2] * m[6];
		m_inv[5] = m[2] * m[3] - m[5] * m[0];
		m_inv[6] = m[3] * m[7] - m[6] * m[4];
		m_inv[7] = m[6] * m[1] - m[0] * m[7];
		m_inv[8] = m[4] * m[0] - m[1] * m[3];
		for (size_t i = 0; i < 9; ++i)
			m_inv[i] /= (2 * tproc.GetVolume(rank) * det);
		Vector3D res = grad_temp;
		res.x *= m_inv[0];
		res.x += grad_temp.y*m_inv[1];
		res.x += grad_temp.z*m_inv[2];
		res.y *= m_inv[3];
		res.y += grad_temp.y*m_inv[4];
		res.y += grad_temp.z*m_inv[5];
		res.z *= m_inv[6];
		res.z += grad_temp.y*m_inv[7];
		res.z += grad_temp.z*m_inv[8];
		return res;
	}
	*/
	Vector3D GetProcCM(Tessellation3D const& tess)
	{
		size_t Ncor = tess.GetPointNo();
		Vector3D res;
		for (size_t i = 0; i < Ncor; ++i)
			res += tess.GetMeshPoint(i);
		res *= 1.0 / static_cast<double>(std::max(static_cast<size_t>(1), Ncor));
		return res;
	}
/*
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
	}*/
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
	/*if (load > 3)
		NewSpeed *= std::min(std::pow(load - 2, 1.5), 0.005 / speed_);
*/
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
	if (tlocal.GetPointNo() == 0)
		RankCM = tproc.GetCellCM(rank);
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
		std::vector<size_t> neigh = tproc.GetNeighbors(rank);
		face_vec faces = tproc.GetCellFaces(rank);
		size_t Nneigh = neigh.size();

		const double neigheps = std::max(0.01,std::min(0.1,RoundSpeed_*0.2));
		for (size_t i = 0; i < Nneigh; ++i)
		{
			if (static_cast<int>(neigh[i]) >= nproc)
				continue;
			const double dist = fastabs(point - tproc.GetMeshPoint(neigh[i]));
			const double mind = neigheps * std::min(MyR, R[neigh[i]]);
			if (dist < mind)
			{
				double speed2 = NewSpeed * 10 * (mind - dist) * std::min(MyR, R[neigh[i]]) / (dist*dist);
				dx += speed2 * (point.x - tproc.GetMeshPoint(neigh[i]).x);
				dy += speed2 * (point.y - tproc.GetMeshPoint(neigh[i]).y);
				dz += speed2 * (point.z - tproc.GetMeshPoint(neigh[i]).z);
			}
			else
			{
				Vector3D otherpoint = RankCMs[neigh[i]];
				//Vector3D otherpoint = tproc.GetMeshPoint(neigh[i]);
				double merit = static_cast<double>(NPerProc[static_cast<size_t>(rank)] - NPerProc[neigh[i]])*std::pow(tproc.GetArea(faces[i]) / (MyR*MyR), 0.4) / IdealPerProc;
				double dr = fastabs(otherpoint - point);
				dx -= NewSpeed * merit*(otherpoint.x - point.x)*MyR / dr;
				dy -= NewSpeed * merit*(otherpoint.y - point.y)*MyR / dr;
				dz -= NewSpeed * merit*(otherpoint.z - point.z)*MyR / dr;
			}
		}
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
