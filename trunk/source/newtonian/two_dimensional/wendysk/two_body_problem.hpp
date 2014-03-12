#ifndef TWO_BODY_PROBLEM_HPP
#define TWO_BODY_PROBLEM_HPP 1

#include "../../../tessellation/geometry.hpp"

using std::pair;

//! \brief Solves the Keplerian problem
class TwoBodyProblem
{
public:

  //! \brief Container for phase space point in 2D
  class GeometricalPoint
  {
  public:

    /*! \brief Class constructor
      \param position_i Position of the point
      \param velocity_i Velocity of the point
     */
    GeometricalPoint(Vector2D const& position_i,
		     Vector2D const& velocity_i);

    virtual ~GeometricalPoint(void);

    //! \brief Position of the point
    const Vector2D position;

    //! \brief Velocity of the point
    const Vector2D velocity;
  };

  //! \brief Container for particle data
  class Particle: public GeometricalPoint
  {
  public:

    /*! \brief Class constructor
      \param position_i Position of the particle
      \param velocity_i Velocity of the particle
      \param mass_i Mass of the particle
     */
    Particle(Vector2D const& position_i,
	     Vector2D const& velocity_i,
	     double mass_i);

    /*! \brief Copy constructor
      \param source Source from which to copy
     */
    Particle(Particle const& source);
    
    //! \brief Mass of the particle
    const double mass;
  };

  /*! \brief Class constructor
    \param pp Initial data of the particles
    \param initial_time Initial time
    \param g Value of the universal gravitation constant
   */
  TwoBodyProblem(pair<Particle, Particle> const& pp,
		 double initial_time=0,
		 double g=1);

  /*! \brief Returns the position and velocities of particles at another time
    \param t Time
    \return Position and velocity of the particles
   */
  std::pair<GeometricalPoint, GeometricalPoint>  timeEvolve(double t) const;  

  /*! \brief Returns the value of the universal gravitation constant
    \return Universal gravitation constant
   */
  double getGravitationConstant(void) const;

  /*! \brief Returns the initial data of the particles
    \return Initial data of the particle
   */
  pair<Particle,Particle> const& getInitialParticles(void) const;

  ~TwoBodyProblem(void);

private:
  const pair<Particle,Particle> pp_;
  const double k_;
  const GeometricalPoint center_of_mass_;
  const GeometricalPoint deviator_;
  const double energy_;
  const double angular_momentum_;
  const Vector2D lrl_;
  const double eccentricity_;
  const double time_scale_;
  const double initial_phase_;
  const double angle_offset_;
};

#endif // TWO_BODY_PROBLEM_HPP
