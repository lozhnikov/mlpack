/**
 * @file hollow_hrect_bound.hpp
 *
 */
#ifndef MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_HPP
#define MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_HPP

#include <mlpack/core.hpp>
#include <mlpack/core/math/range.hpp>
#include <mlpack/core/metrics/lmetric.hpp>
#include "hrectbound.hpp"
#include "bound_traits.hpp"

namespace mlpack {
namespace bound {

/**
 * @tparam MetricType Type of metric to use; must be of type LMetric.
 * @tparam ElemType Element type (double/float/int/etc.).
 */
template<typename TMetricType = metric::LMetric<2, true>,
         typename ElemType = double>
class HollowHRectBound
{
 public:
 typedef TMetricType MetricType;
  /**
   * Empty constructor; creates a bound of dimensionality 0.
   */
  HollowHRectBound();

  /**
   * Initializes to specified dimensionality with each dimension the empty
   * set.
   */
  HollowHRectBound(const size_t dimension);

  //! Copy constructor; necessary to prevent memory leaks.
  HollowHRectBound(const HollowHRectBound& other);

  inline HollowHRectBound<TMetricType, ElemType>& operator=(
    const HollowHRectBound<TMetricType, ElemType>& other);

  //! Move constructor: take possession of another bound's information.
  HollowHRectBound(HollowHRectBound&& other);

  /**
   * Resets all dimensions to the empty set (so that this bound contains
   * nothing).
   */
  void Clear();

  //! Gets the dimensionality.
  size_t Dim() const { return outerRect.Dim(); }

  //! Get the range for a particular dimension.  No bounds checking.  Be
  //! careful: this may make MinWidth() invalid.
  math::RangeType<ElemType>& operator[](const size_t i) { return outerRect[i]; }
  //! Modify the range for a particular dimension.  No bounds checking.
  const math::RangeType<ElemType>& operator[](const size_t i) const
  { return outerRect[i]; }

  //! Get the minimum width of the bound.
  ElemType MinWidth() const { return outerRect.MinWidth(); }
  //! Modify the minimum width of the bound.
  ElemType& MinWidth() { return outerRect.MinWidth(); }

  const HRectBound<MetricType, ElemType>& OuterRect() const
  { return outerRect; }

  std::list<HRectBound<MetricType, ElemType>>& Hollows() { return hollows; }

  const std::list<HRectBound<MetricType, ElemType>>& Hollows() const
  { return hollows; }

  const MetricType& Metric() const { return *metric; }


  /**
   * Calculates the center of the range, placing it into the given vector.
   *
   * @param center Vector which the center will be written to.
   */
  void Center(arma::Col<ElemType>& center) const { outerRect.Center(center); };

  /**
   * Calculates minimum bound-to-point distance.
   *
   * @param point Point to which the minimum distance is requested.
   */
  template<typename VecType>
  ElemType MinDistance(const VecType& point,
                       typename boost::enable_if<IsVector<VecType>>* = 0) const;

  /**
   * Calculates minimum bound-to-bound distance.
   *
   * @param other Bound to which the minimum distance is requested.
   */
  ElemType MinDistance(const HollowHRectBound& other) const;

  /**
   * Calculates maximum bound-to-point squared distance.
   *
   * @param point Point to which the maximum distance is requested.
   */
  template<typename VecType>
  ElemType MaxDistance(const VecType& point,
                       typename boost::enable_if<IsVector<VecType>>* = 0) const;

  /**
   * Computes maximum distance.
   *
   * @param other Bound to which the maximum distance is requested.
   */
  ElemType MaxDistance(const HollowHRectBound& other) const;

  /**
   * Calculates minimum and maximum bound-to-bound distance.
   *
   * @param other Bound to which the minimum and maximum distances are
   *     requested.
   */
  math::RangeType<ElemType> RangeDistance(const HollowHRectBound& other) const;

  /**
   * Calculates minimum and maximum bound-to-point distance.
   *
   * @param point Point to which the minimum and maximum distances are
   *     requested.
   */
  template<typename VecType>
  math::RangeType<ElemType> RangeDistance(
      const VecType& point,
      typename boost::enable_if<IsVector<VecType>>* = 0) const;

  /**
   * Expands this region to include new points.
   *
   * @tparam MatType Type of matrix; could be Mat, SpMat, a subview, or just a
   *   vector.
   * @param data Data points to expand this region to include.
   */
  template<typename MatType>
  HollowHRectBound& operator|=(const MatType& data);

  /**
   * Expands this region to encompass another bound.
   */
  HollowHRectBound& operator|=(const HollowHRectBound& other);

  /**
   * Determines if a point is within this bound.
   */
  template<typename VecType>
  bool Contains(const VecType& point) const;

  /**
   * Determines if this bound partially contains a bound.
   */
  bool Contains(const HollowHRectBound& bound) const;

  /**
   * Returns the diameter of the hyperrectangle (that is, the longest diagonal).
   */
  ElemType Diameter() const { return outerRect.Diameter(); };

  void AddHollow(const HRectBound<MetricType, ElemType>& hollow);

  /**
   * Serialize the bound object.
   */
  template<typename Archive>
  void Serialize(Archive& ar, const unsigned int version);

 private:
  HRectBound<MetricType, ElemType> outerRect;
  std::list<HRectBound<MetricType, ElemType>> hollows;
  static constexpr size_t maxNumHollows = 10;
  //! The metric used in this bound.
  MetricType* metric;

  /**
   * To know whether this object allocated memory to the metric member
   * variable. This will be true except in the copy constructor and the
   * overloaded assignment operator. We need this to know whether we should
   * delete the metric member variable in the destructor.
   */
  bool ownsMetric;
};

// A specialization of BoundTraits for this class.
template<typename MetricType, typename ElemType>
struct BoundTraits<HollowHRectBound<MetricType, ElemType>>
{
  //! These bounds are always tight for each dimension.
  const static bool HasTightBounds = true;
};

namespace meta {
template<typename BoundType>
struct IsHollowHRectBound
{
  static const bool value = false;
};

template<typename MetricType, typename ElemType>
struct IsHollowHRectBound<HollowHRectBound<MetricType, ElemType>>
{
  static const bool value = true;
};

} // namespace meta

} // namespace bound
} // namespace mlpack

#include "hollow_hrect_bound_impl.hpp"

#endif // MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_HPP

