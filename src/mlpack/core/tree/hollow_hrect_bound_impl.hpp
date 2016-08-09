/**
 * @file hollow_hrect_bound_impl.hpp
 *
 */
#ifndef MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_IMPL_HPP
#define MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_IMPL_HPP

#include <math.h>

// In case it has not been included yet.
#include "hollow_hrect_bound.hpp"

namespace mlpack {
namespace bound {

/**
 * Empty constructor.
 */
template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>::HollowHRectBound():
    metric(new MetricType()),
    ownsMetric(true)

{ /* Nothing to do. */ }

/**
 * Initializes to specified dimensionality with each dimension the empty
 * set.
 */
template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>::
HollowHRectBound(const size_t dimension) :
    outerRect(dimension),
    metric(new MetricType()),
    ownsMetric(true)
{ /* Nothing to do. */ }

/**
 * Copy constructor necessary to prevent memory leaks.
 */
template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>::HollowHRectBound(
    const HollowHRectBound<TMetricType, ElemType>& other) :
    outerRect(other.outerRect),
    hollows(other.hollows),
    metric(other.metric),
    ownsMetric(false)
{

}

template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>&
HollowHRectBound<TMetricType, ElemType>::operator=(
    const HollowHRectBound<TMetricType, ElemType>& other)
{
  outerRect = other.outerRect;
  hollows = other.hollows;
  metric = other.metric;
  ownsMetric = false;
  return *this;
}

/**
 * Move constructor: take possession of another bound's information.
 */
template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>::HollowHRectBound(
    HollowHRectBound<TMetricType, ElemType>&& other) :
    outerRect(std::move(other.outerRect)),
    hollows(std::move(other.hollows)),
    metric(other.metric),
    ownsMetric(other.ownsMetric)

{
  other.metric = NULL;
  other.ownsMetric = false;
}

/**
 * Resets all dimensions to the empty set.
 */
template<typename TMetricType, typename ElemType>
inline void HollowHRectBound<TMetricType, ElemType>::Clear()
{
  outerRect.Clear();
  hollows.clear();
}

/**
 * Calculates minimum bound-to-point squared distance.
 */
template<typename TMetricType, typename ElemType>
template<typename VecType>
inline ElemType HollowHRectBound<TMetricType, ElemType>::
MinDistance(const VecType& point,
            typename boost::enable_if<IsVector<VecType>>* /* junk */) const
{
  ElemType dist = outerRect.MinDistance(point);

  if (dist == 0 && hollows.size() > 0)
  {
    typedef typename std::list<HRectBound<MetricType, ElemType>>::const_iterator iterator;

    for (iterator it = hollows.begin(); it != hollows.end(); it++)
    {
      ElemType hollowDist = math::ClampNonNegative(point[0] - (*it)[0].Lo());
      for (size_t k = 0; k < outerRect.Dim(); k++)
      {
        const ElemType loDist = math::ClampNonNegative(point[k] - (*it)[k].Lo());
        const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() - point[k]);

        hollowDist = 0.5 * (hollowDist + loDist -
            std::fabs(hollowDist - loDist));
        hollowDist = 0.5 * (hollowDist + hiDist -
            std::fabs(hollowDist - hiDist));
      }
      dist = 0.5 * (dist + hollowDist + std::fabs(dist - hollowDist));
    }
  }

  return dist;
}

/**
 * Calculates minimum bound-to-bound squared distance.
 */
template<typename TMetricType, typename ElemType>
ElemType HollowHRectBound<TMetricType, ElemType>::MinDistance(
    const HollowHRectBound& other) const
{
  ElemType dist = outerRect.MinDistance(other.outerRect);

  if (dist == 0 && hollows.size() > 0)
  {
    typedef typename std::list<HRectBound<MetricType, ElemType>>::const_iterator iterator;
    bool firstIncludedInSecond = true;
    bool secondIncludedInFirst = true;

    for (size_t k = 0; k < outerRect.Dim(); k++)
    {
      if (outerRect[k].Lo() < other.outerRect[k].Lo())
        firstIncludedInSecond = false;
      else if (outerRect[k].Lo() > other.outerRect[k].Lo())
        secondIncludedInFirst = false;
    }

    if (!firstIncludedInSecond && !secondIncludedInFirst)
      return 0;
    else if (firstIncludedInSecond)
    {
      for (iterator it = other.hollows.begin(); it != other.hollows.end(); it++)
      {
        ElemType hollowDist = math::ClampNonNegative(outerRect[0].Lo() -
              (*it)[0].Lo());
        for (size_t k = 0; k < outerRect.Dim(); k++)
        {
          const ElemType loDist = math::ClampNonNegative(outerRect[k].Lo() -
              (*it)[k].Lo());
          const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() -
              outerRect[k].Hi());

          hollowDist = 0.5 * (hollowDist + loDist -
              std::fabs(hollowDist - loDist));
          hollowDist = 0.5 * (hollowDist + hiDist -
              std::fabs(hollowDist - hiDist));
        }
        dist = 0.5 * (dist + hollowDist + std::fabs(dist - hollowDist));
      }
    }
    else if (secondIncludedInFirst)
    {
      for (iterator it = hollows.begin(); it != hollows.end(); it++)
      {
        ElemType hollowDist = math::ClampNonNegative(other.outerRect[0].Lo() -
              (*it)[0].Lo());
        for (size_t k = 0; k < other.outerRect.Dim(); k++)
        {
          const ElemType loDist = math::ClampNonNegative(other.outerRect[k].Lo() -
              (*it)[k].Lo());
          const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() -
              other.outerRect[k].Hi());

          hollowDist = 0.5 * (hollowDist + loDist -
              std::fabs(hollowDist - loDist));
          hollowDist = 0.5 * (hollowDist + hiDist -
              std::fabs(hollowDist - hiDist));
        }
        dist = 0.5 * (dist + hollowDist + std::fabs(dist - hollowDist));
      }
    }
  }

  return dist;
}

/**
 * Calculates maximum bound-to-point squared distance.
 */
template<typename TMetricType, typename ElemType>
template<typename VecType>
inline ElemType HollowHRectBound<TMetricType, ElemType>::
MaxDistance(const VecType& point,
            typename boost::enable_if<IsVector<VecType> >* /* junk */) const
{
  return outerRect.MaxDistance(point);
}

/**
 * Computes maximum distance.
 */
template<typename TMetricType, typename ElemType>
inline ElemType HollowHRectBound<TMetricType, ElemType>::
MaxDistance(const HollowHRectBound& other) const
{
  return outerRect.MaxDistance(other.outerRect);
}

/**
 * Calculates minimum and maximum bound-to-bound squared distance.
 */
template<typename TMetricType, typename ElemType>
inline math::RangeType<ElemType>
HollowHRectBound<TMetricType, ElemType>::RangeDistance(
    const HollowHRectBound& other) const
{
  math::RangeType<ElemType> range = outerRect.RangeDistance(other.outerRect);

  if (range.Lo() == 0 && hollows.size() > 0)
  {
    typedef typename std::list<HRectBound<MetricType, ElemType>>::const_iterator iterator;
    bool firstIncludedInSecond = true;
    bool secondIncludedInFirst = true;

    for (size_t k = 0; k < outerRect.Dim(); k++)
    {
      if (outerRect[k].Lo() < other.outerRect[k].Lo())
        firstIncludedInSecond = false;
      else if (outerRect[k].Lo() > other.outerRect[k].Lo())
        secondIncludedInFirst = false;
    }

    if (!firstIncludedInSecond && !secondIncludedInFirst)
      return range;
    else if (firstIncludedInSecond)
    {
      for (iterator it = other.hollows.begin(); it != other.hollows.end(); it++)
      {
        ElemType hollowDist = math::ClampNonNegative(outerRect[0].Lo() -
              (*it)[0].Lo());
        for (size_t k = 0; k < outerRect.Dim(); k++)
        {
          const ElemType loDist = math::ClampNonNegative(outerRect[k].Lo() -
              (*it)[k].Lo());
          const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() -
              outerRect[k].Hi());

          hollowDist = 0.5 * (hollowDist + loDist -
              std::fabs(hollowDist - loDist));
          hollowDist = 0.5 * (hollowDist + hiDist -
              std::fabs(hollowDist - hiDist));
        }
        range.Lo() = 0.5 * (range.Lo() + hollowDist +
            std::fabs(range.Lo() - hollowDist));
      }
    }
    else if (secondIncludedInFirst)
    {
      for (iterator it = hollows.begin(); it != hollows.end(); it++)
      {
        ElemType hollowDist = math::ClampNonNegative(other.outerRect[0].Lo() -
              (*it)[0].Lo());
        for (size_t k = 0; k < other.outerRect.Dim(); k++)
        {
          const ElemType loDist = math::ClampNonNegative(other.outerRect[k].Lo() -
              (*it)[k].Lo());
          const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() -
              other.outerRect[k].Hi());

          hollowDist = 0.5 * (hollowDist + loDist -
              std::fabs(hollowDist - loDist));
          hollowDist = 0.5 * (hollowDist + hiDist -
              std::fabs(hollowDist - hiDist));
        }
        range.Lo() = 0.5 * (range.Lo() + hollowDist +
            std::fabs(range.Lo() - hollowDist));
      }
    }
  }
  return range;
}

/**
 * Calculates minimum and maximum bound-to-point squared distance.
 */
template<typename TMetricType, typename ElemType>
template<typename VecType>
inline math::RangeType<ElemType>
HollowHRectBound<TMetricType, ElemType>::RangeDistance(
    const VecType& point,
    typename boost::enable_if<IsVector<VecType>>* /* junk */) const
{
  math::RangeType<ElemType> range = outerRect.RangeDistance(point);

  if (range.Lo() == 0 && hollows.size() > 0)
  {
    typedef typename std::list<HRectBound<MetricType, ElemType>>::const_iterator iterator;

    for (iterator it = hollows.begin(); it != hollows.end(); it++)
    {
      ElemType hollowDist = math::ClampNonNegative(point[0] - (*it)[0].Lo());;
      for (size_t k = 0; k < outerRect.Dim(); k++)
      {
        const ElemType loDist = math::ClampNonNegative(point[k] - (*it)[k].Lo());
        const ElemType hiDist = math::ClampNonNegative((*it)[k].Hi() - point[k]);

        hollowDist = 0.5 * (hollowDist + loDist -
            std::fabs(hollowDist - loDist));
        hollowDist = 0.5 * (hollowDist + hiDist -
            std::fabs(hollowDist - hiDist));
      }
      range.Lo() = 0.5 * (range.Lo() + hollowDist +
          std::fabs(range.Lo() - hollowDist));
    }
  }

  return range;
}

/**
 * Expands this region to include a new point.
 */
template<typename TMetricType, typename ElemType>
template<typename MatType>
inline HollowHRectBound<TMetricType, ElemType>&
HollowHRectBound<TMetricType, ElemType>::operator|=(
    const MatType& data)
{
  outerRect |= data;
  hollows.clear();
  return *this;
}

/**
 * Expands this region to encompass another bound.
 */
template<typename TMetricType, typename ElemType>
inline HollowHRectBound<TMetricType, ElemType>&
HollowHRectBound<TMetricType, ElemType>::operator|=(
    const HollowHRectBound& other)
{
  outerRect |= other;
  hollows.clear();

  return *this;
}

/**
 * Determines if a point is within this bound.
 */
template<typename TMetricType, typename ElemType>
template<typename VecType>
inline bool HollowHRectBound<TMetricType, ElemType>::Contains(
    const VecType& point) const
{
  if (!outerRect.Contains(point))
    return false;

  typedef typename std::list<HRectBound<MetricType, ElemType>>::const_iterator iterator;

  for (iterator it = hollows.begin(); it != hollows.end(); it++)
  {
    bool success = false;
    for (size_t k = 0; k < outerRect.Dim(); k++)
      if ( (*it)[k].Lo() >= point[k] || (*it)[k].Hi() <= point[k])
        success = true;
    if (!success)
      return false;
  }
  return true;
}

template<typename TMetricType, typename ElemType>
inline void HollowHRectBound<TMetricType, ElemType>::AddHollow(
    const HRectBound<MetricType, ElemType>& hollow)
{
  if (hollows.size() < maxNumHollows)
  {
    hollows.push_back(hollow);
    return;
  }

/*  ElemType volume = hollow.Volume();

  typedef typename std::list<HRectBound<MetricType, ElemType>>::iterator iterator;
  iterator minVolumeHollow;
  ElemType minVolume = std::numeric_limits<ElemType>::max();

  for (iterator it = hollows.begin(); it != hollows.end(); it++)
  {
    const ElemType vol = (*it).Volume();

    if (vol < minVolume)
    {
      minVolume = vol;
      minVolumeHollow = it;
    }
  }

  if (minVolume < volume)
    *minVolumeHollow = hollow;
*/

  typedef typename std::list<HRectBound<MetricType, ElemType>>::iterator iterator;
  iterator minVolumeHollow;
  HRectBound<MetricType, ElemType> bound(Dim());

  for (iterator it = hollows.begin(); it != hollows.end(); it++)
    bound |= *it;

  ElemType volume = bound.Volume();
  ElemType maxVolume = 0;

  for (iterator it = hollows.begin(); it != hollows.end(); it++)
  {
    bound = hollow;

    for (iterator it2 = hollows.begin(); it2 != hollows.end(); it2++)
    {
      if (it == it2)
        continue;
      bound |= *it2;
    }
    ElemType vol = bound.Volume();
    if (vol > maxVolume)
    {
      maxVolume = vol;
      minVolumeHollow = it;
    }    
  }

  if (maxVolume > volume)
    *minVolumeHollow = hollow;
}

//! Serialize the bound object.
template<typename TMetricType, typename ElemType>
template<typename Archive>
void HollowHRectBound<TMetricType, ElemType>::Serialize(
    Archive& ar, const unsigned int /* version */)
{
  ar & data::CreateNVP(outerRect, "outerRect");
//  ar & data::CreateNVP(hollows, "hollows");
}

} // namespace bound
} // namespace mlpack

#endif // MLPACK_CORE_TREE_HOLLOW_HRECT_BOUND_IMPL_HPP

