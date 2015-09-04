/*! \file simple_cell_updater.hpp
  \author Almog Yalinewich
  \brief Simple cell updater
*/

#ifndef SIMPLE_CELL_UPDATER_HPP
#define SIMPLE_CELL_UPDATER_HPP 1

#include "cell_updater_2d.hpp"
#include "../../misc/utils.hpp"

using std::vector;

//! \brief Simple cell updater
class SimpleCellUpdater: public CellUpdater
{
public:

  class Condition
  {
  public:

    virtual ~Condition(void);

    virtual bool operator()
    (const Tessellation& tess,
     const PhysicalGeometry& pg,
     const EquationOfState& eos,
     const vector<Extensive>& extensives,
     const vector<ComputationalCell>& cells,
     const CacheData& cd,
     const size_t index) const = 0;
  };

  class Action
  {
  public:

    virtual ~Action(void);

    virtual ComputationalCell operator()
    (const Tessellation& tess,
     const PhysicalGeometry& pg,
     const EquationOfState& eos,
     const vector<Extensive>& extensives,
     const vector<ComputationalCell>& cells,
     const CacheData& cd,
     const size_t index) const = 0;
  };

  SimpleCellUpdater
  (const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence =
   vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >());

  vector<ComputationalCell> operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const;

  ~SimpleCellUpdater(void);

private:
  const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence_;
};

class HasSticker: public SimpleCellUpdater::Condition
{
public:

  HasSticker(const string& sticker_name);

  bool operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& cells,
   const CacheData& cd,
   const size_t index) const;

private:
  const string sticker_name_;
};

class SkipUpdate: public SimpleCellUpdater::Action
{
public:

  SkipUpdate(void);

  ComputationalCell operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& cells,
   const CacheData& cd,
   const size_t index) const;
};

#endif // SIMPLE_CELL_UPDATER_HPP
