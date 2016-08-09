/**
 * @file tree_test.cpp
 *
 * Tests for tree-building methods.
 */
#include <mlpack/core.hpp>
#include <mlpack/core/tree/bounds.hpp>
#include <mlpack/core/metrics/lmetric.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <mlpack/core/tree/vantage_point_tree.hpp>

#include <boost/test/unit_test.hpp>
#include "test_tools.hpp"

using namespace mlpack;
using namespace mlpack::math;
using namespace mlpack::tree;
using namespace mlpack::neighbor;
using namespace mlpack::metric;
using namespace mlpack::bound;

BOOST_AUTO_TEST_SUITE(VantagePointTreeTest);

BOOST_AUTO_TEST_CASE(VPTreeTraitsTest)
{
  typedef VPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;

  bool b = TreeTraits<TreeType>::HasOverlappingChildren;
  BOOST_REQUIRE_EQUAL(b, true);
  b = TreeTraits<TreeType>::FirstPointIsCentroid;
  BOOST_REQUIRE_EQUAL(b, false);
  b = TreeTraits<TreeType>::HasSelfChildren;
  BOOST_REQUIRE_EQUAL(b, false);
  b = TreeTraits<TreeType>::RearrangesDataset;
  BOOST_REQUIRE_EQUAL(b, true);
  b = TreeTraits<TreeType>::BinaryTree;
  BOOST_REQUIRE_EQUAL(b, false);
}

BOOST_AUTO_TEST_CASE(HollowBallBoundTest)
{
  HollowBallBound<EuclideanDistance> b(2, 4, arma::vec("1.0 2.0 3.0 4.0 5.0"));

  BOOST_REQUIRE_EQUAL(b.Contains(arma::vec("1.0 2.0 3.0 7.0 5.0")), true);

  BOOST_REQUIRE_EQUAL(b.Contains(arma::vec("1.0 2.0 3.0 9.0 5.0")), false);

  BOOST_REQUIRE_EQUAL(b.Contains(arma::vec("1.0 2.0 3.0 5.0 5.0")), false);

  HollowBallBound<EuclideanDistance> b2(0.5, 1,
      arma::vec("1.0 2.0 3.0 7.0 5.0"));
  BOOST_REQUIRE_EQUAL(b.Contains(b2), true);

  b2 = HollowBallBound<EuclideanDistance>(2.5, 3.5,
      arma::vec("1.0 2.0 3.0 4.5 5.0"));
  BOOST_REQUIRE_EQUAL(b.Contains(b2), true);

  b2 = HollowBallBound<EuclideanDistance>(2.0, 3.5,
      arma::vec("1.0 2.0 3.0 4.5 5.0"));
  BOOST_REQUIRE_EQUAL(b.Contains(b2), false);

  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("1.0 2.0 8.0 4.0 5.0")), 1.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("1.0 2.0 4.0 4.0 5.0")), 1.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("1.0 2.0 3.0 4.0 5.0")), 2.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("1.0 2.0 5.0 4.0 5.0")), 0.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("5.0 2.0 3.0 4.0 5.0")), 0.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MinDistance(arma::vec("3.0 2.0 3.0 4.0 5.0")), 0.0,
      1e-5);

  BOOST_REQUIRE_CLOSE(b.MaxDistance(arma::vec("1.0 2.0 4.0 4.0 5.0")), 5.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MaxDistance(arma::vec("1.0 2.0 8.0 4.0 5.0")), 9.0,
      1e-5);
  BOOST_REQUIRE_CLOSE(b.MaxDistance(arma::vec("1.0 2.0 3.0 4.0 5.0")), 4.0,
      1e-5);

  b2 = HollowBallBound<EuclideanDistance>(3, 4,
      arma::vec("1.0 2.0 3.0 5.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MinDistance(b2), 0.0, 1e-5);

  b2 = HollowBallBound<EuclideanDistance>(1, 2,
      arma::vec("1.0 2.0 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MinDistance(b2), 0.0, 1e-5);

  b2 = HollowBallBound<EuclideanDistance>(0.5, 1.0,
      arma::vec("1.0 2.5 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MinDistance(b2), 0.5, 1e-5);

  b2 = HollowBallBound<EuclideanDistance>(0.5, 1.0,
      arma::vec("1.0 8.0 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MinDistance(b2), 1.0, 1e-5);

  b2 = HollowBallBound<EuclideanDistance>(0.5, 2.0,
      arma::vec("1.0 8.0 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MinDistance(b2), 0.0, 1e-5);

  b2 = HollowBallBound<EuclideanDistance>(0.5, 2.0,
      arma::vec("1.0 8.0 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MaxDistance(b2), 12.0, 1e-5);
  
  b2 = HollowBallBound<EuclideanDistance>(0.5, 2.0,
      arma::vec("1.0 3.0 3.0 4.0 5.0"));
  BOOST_REQUIRE_CLOSE(b.MaxDistance(b2), 7.0, 1e-5);

  HollowBallBound<EuclideanDistance> b1 = b;
  b2 = HollowBallBound<EuclideanDistance>(1.0, 2.0,
      arma::vec("1.0 2.5 3.0 4.0 5.0"));

  b1 |= b2;
  BOOST_REQUIRE_CLOSE(b1.InnerRadius(), 0.5, 1e-5);
  
  b1 = b;
  b2 = HollowBallBound<EuclideanDistance>(0.5, 2.0,
      arma::vec("1.0 3.0 3.0 4.0 5.0"));
  b1 |= b2;
  BOOST_REQUIRE_CLOSE(b1.InnerRadius(), 0.0, 1e-5);

  b1 = b;
  b2 = HollowBallBound<EuclideanDistance>(0.5, 4.0,
      arma::vec("1.0 3.0 3.0 4.0 5.0"));
  b1 |= b2;
  BOOST_REQUIRE_CLOSE(b1.OuterRadius(), 5.0, 1e-5);
}

template<typename TreeType>
void CheckBound(TreeType& tree)
{
  typedef typename TreeType::ElemType ElemType;
  if (tree.IsLeaf())
  {
    // Ensure that the bound contains all descendant points.
    for (size_t i = 0; i < tree.NumPoints(); i++)
    {
      ElemType dist = tree.Bound().Metric().Evaluate(tree.Bound().Center(),
        tree.Dataset().col(tree.Point(i)));

      BOOST_REQUIRE_LE(tree.Bound().InnerRadius(), dist  *
          (1.0 + 10.0 * std::numeric_limits<ElemType>::epsilon()));

      BOOST_REQUIRE_LE(dist, tree.Bound().OuterRadius() *
          (1.0 + 10.0 * std::numeric_limits<ElemType>::epsilon()));
    }
  }
  else
  {
    TreeType* central = tree.Central();
    BOOST_REQUIRE_EQUAL(central->NumPoints(), 1);
    BOOST_REQUIRE_EQUAL(true,
        central->Bound().Contains(tree.Dataset().col(central->Point(0))));
    BOOST_REQUIRE_EQUAL(central->Bound().InnerRadius(), 0.0);
    BOOST_REQUIRE_EQUAL(central->Bound().OuterRadius(), 0.0);

    // Ensure that the bound contains all descendant points.
    for (size_t i = 0; i < tree.NumDescendants(); i++)
    {
      ElemType dist = tree.Bound().Metric().Evaluate(tree.Bound().Center(),
        tree.Dataset().col(tree.Descendant(i)));

      BOOST_REQUIRE_LE(tree.Bound().InnerRadius(), dist  *
          (1.0 + 10.0 * std::numeric_limits<ElemType>::epsilon()));

      BOOST_REQUIRE_LE(dist, tree.Bound().OuterRadius() *
          (1.0 + 10.0 * std::numeric_limits<ElemType>::epsilon()));
    }

    CheckBound(*tree.Inner());
    CheckBound(*tree.Outer());
  }
}

BOOST_AUTO_TEST_CASE(VPTreeBoundTest)
{
  typedef VPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;

  arma::mat dataset(8, 1000);
  dataset.randu();

  TreeType tree(dataset);
  CheckBound(tree);
}

template<typename TreeType>
void CheckSplit(TreeType& tree)
{
  if(tree.IsLeaf())
    return;

  typename TreeType::ElemType maxDist = 0;

  size_t pointsEnd = tree.Inner()->Begin() + tree.Inner()->Count();
  for (size_t i = tree.Inner()->Begin(); i < pointsEnd; i++)
  {
    typename TreeType::ElemType dist =
        tree.Bound().Metric().Evaluate(tree.Dataset().col(i),
            tree.Dataset().col(tree.Central()->Begin()));

    if (dist > maxDist)
      maxDist = dist;
  }

  pointsEnd = tree.Outer()->Begin() + tree.Outer()->Count();
  for (size_t i = tree.Outer()->Begin(); i < pointsEnd; i++)
  {
    typename TreeType::ElemType dist =
        tree.Bound().Metric().Evaluate(tree.Dataset().col(i),
            tree.Dataset().col(tree.Central()->Begin()));
    BOOST_REQUIRE_LE(maxDist, dist);
  }

  for (size_t k = 0; k < tree.Bound().Dim(); k++)
  {
    BOOST_REQUIRE_EQUAL(tree.Inner()->Bound().Center()[k],
        tree.Dataset().col(tree.Central()->Point(0))[k]);
    BOOST_REQUIRE_EQUAL(tree.Outer()->Bound().Center()[k],
        tree.Dataset().col(tree.Central()->Point(0))[k]);
    BOOST_REQUIRE_EQUAL(tree.Central()->Bound().Center()[k],
        tree.Dataset().col(tree.Central()->Point(0))[k]);
  }

  CheckSplit(*tree.Inner());
  CheckSplit(*tree.Outer());
}

BOOST_AUTO_TEST_CASE(VPTreeSplitTest)
{
  typedef VPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;

  arma::mat dataset(8, 1000);
  dataset.randu();

  TreeType tree(dataset);
  CheckSplit(tree);
}

BOOST_AUTO_TEST_CASE(VPTreeTest)
{
  typedef VPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;

  size_t maxRuns = 10; // Ten total tests.
  size_t pointIncrements = 1000; // Range is from 2000 points to 11000.

  // We use the default leaf size of 20.
  for (size_t run = 0; run < maxRuns; run++)
  {
    size_t dimensions = run + 2;
    size_t maxPoints = (run + 1) * pointIncrements;

    size_t size = maxPoints;
    arma::mat dataset = arma::mat(dimensions, size);
    arma::mat datacopy; // Used to test mappings.

    // Mappings for post-sort verification of data.
    std::vector<size_t> newToOld;
    std::vector<size_t> oldToNew;

    // Generate data.
    dataset.randu();

    // Build the tree itself.
    TreeType root(dataset, newToOld, oldToNew);
    const arma::mat& treeset = root.Dataset();

    // Ensure the size of the tree is correct.
    BOOST_REQUIRE_EQUAL(root.NumDescendants(), size);

    // Check the forward and backward mappings for correctness.
    for(size_t i = 0; i < size; i++)
    {
      for(size_t j = 0; j < dimensions; j++)
      {
        BOOST_REQUIRE_EQUAL(treeset(j, i), dataset(j, newToOld[i]));
        BOOST_REQUIRE_EQUAL(treeset(j, oldToNew[i]), dataset(j, i));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(SingleTreeTraverserTest)
{
  arma::mat dataset;
  dataset.randu(8, 1000); // 1000 points in 8 dimensions.
  arma::Mat<size_t> neighbors1;
  arma::mat distances1;
  arma::Mat<size_t> neighbors2;
  arma::mat distances2;

  // Nearest neighbor search with the VP tree.
  NeighborSearch<NearestNeighborSort, metric::LMetric<2, true>, arma::mat,
      VPTree> knn1(dataset, false, true);

  knn1.Search(5, neighbors1, distances1);

  // Nearest neighbor search the naive way.
  KNN knn2(dataset, true, true);

  knn2.Search(5, neighbors2, distances2);

  for (size_t i = 0; i < neighbors1.size(); i++)
  {
    BOOST_REQUIRE_EQUAL(neighbors1[i], neighbors2[i]);
    BOOST_REQUIRE_EQUAL(distances1[i], distances2[i]);
  }
}

BOOST_AUTO_TEST_CASE(DualTreeTraverserTest)
{
  arma::mat dataset;
  dataset.randu(8, 1000); // 1000 points in 8 dimensions.
  arma::Mat<size_t> neighbors1;
  arma::mat distances1;
  arma::Mat<size_t> neighbors2;
  arma::mat distances2;

  // Nearest neighbor search with the VP tree.
  NeighborSearch<NearestNeighborSort, metric::LMetric<2, true>, arma::mat,
      VPTree> knn1(dataset, false, false);

  knn1.Search(5, neighbors1, distances1);

  // Nearest neighbor search the naive way.
  KNN knn2(dataset, true, true);

  knn2.Search(5, neighbors2, distances2);

  for (size_t i = 0; i < neighbors1.size(); i++)
  {
    BOOST_REQUIRE_EQUAL(neighbors1[i], neighbors2[i]);
    BOOST_REQUIRE_EQUAL(distances1[i], distances2[i]);
  }
}

BOOST_AUTO_TEST_CASE(HRectVPTreeTest)
{
  typedef HRectVPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;

  size_t maxRuns = 10; // Ten total tests.
  size_t pointIncrements = 1000; // Range is from 2000 points to 11000.

  // We use the default leaf size of 20.
  for (size_t run = 0; run < maxRuns; run++)
  {
    size_t dimensions = run + 2;
    size_t maxPoints = (run + 1) * pointIncrements;

    size_t size = maxPoints;
    arma::mat dataset = arma::mat(dimensions, size);
    arma::mat datacopy; // Used to test mappings.

    // Mappings for post-sort verification of data.
    std::vector<size_t> newToOld;
    std::vector<size_t> oldToNew;

    // Generate data.
    dataset.randu();

    // Build the tree itself.
    TreeType root(dataset, newToOld, oldToNew);
    const arma::mat& treeset = root.Dataset();

    // Ensure the size of the tree is correct.
    BOOST_REQUIRE_EQUAL(root.NumDescendants(), size);

    // Check the forward and backward mappings for correctness.
    for(size_t i = 0; i < size; i++)
    {
      for(size_t j = 0; j < dimensions; j++)
      {
        BOOST_REQUIRE_EQUAL(treeset(j, i), dataset(j, newToOld[i]));
        BOOST_REQUIRE_EQUAL(treeset(j, oldToNew[i]), dataset(j, i));
      }
    }
  }
}

template<typename TreeType, typename MetricType>
void CheckDistance(TreeType& tree, TreeType* node = NULL)
{
  typedef typename TreeType::ElemType ElemType;
  if (node == NULL)
  {
    node = &tree;

    while (node->Parent() != NULL)
      node = node->Parent();

    CheckDistance<TreeType, MetricType>(tree, node);

    for (size_t j = 0; j < tree.Dataset().n_cols; j++)
    {
      const arma::Col<ElemType>& point = tree.  Dataset().col(j);
      ElemType maxDist = 0;
      ElemType minDist = std::numeric_limits<ElemType>::max();
      for (size_t i = 0; i < tree.NumDescendants(); i++)
      {
        ElemType dist = MetricType::Evaluate(
            tree.Dataset().col(tree.Descendant(i)),
            tree.Dataset().col(j));

        if (dist > maxDist)
          maxDist = dist;
        if (dist < minDist)
          minDist = dist;
      }

      BOOST_REQUIRE_LE(tree.Bound().MinDistance(point), minDist *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
      BOOST_REQUIRE_LE(maxDist, tree.Bound().MaxDistance(point) *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));

      math::RangeType<ElemType> r = tree.Bound().RangeDistance(point);

      BOOST_REQUIRE_LE(r.Lo(), minDist *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
      BOOST_REQUIRE_LE(maxDist, r.Hi() *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
    }
      
    if (!tree.IsLeaf())
    {
      CheckDistance<TreeType, MetricType>(*tree.Central());
      CheckDistance<TreeType, MetricType>(*tree.Inner());
      CheckDistance<TreeType, MetricType>(*tree.Outer());
    }
  }
  else
  {
    if (&tree != node)
    {
      ElemType maxDist = 0;
      ElemType minDist = std::numeric_limits<ElemType>::max();
      for (size_t i = 0; i < tree.NumDescendants(); i++)
        for (size_t j = 0; j < node->NumDescendants(); j++)
        {
          ElemType dist = MetricType::Evaluate(
              tree.Dataset().col(tree.Descendant(i)),
              node->Dataset().col(node->Descendant(j)));

          if (dist > maxDist)
            maxDist = dist;
          if (dist < minDist)
            minDist = dist;
        }

      BOOST_REQUIRE_LE(tree.Bound().MinDistance(node->Bound()), minDist *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
      BOOST_REQUIRE_LE(maxDist, tree.Bound().MaxDistance(node->Bound()) *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));

      math::RangeType<ElemType> r = tree.Bound().RangeDistance(node->Bound());

      BOOST_REQUIRE_LE(r.Lo(), minDist *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
      BOOST_REQUIRE_LE(maxDist, r.Hi() *
          (1.0 + 10 * std::numeric_limits<ElemType>::epsilon()));
    }
    if (!node->IsLeaf())
    {
      CheckDistance<TreeType, MetricType>(tree, node->Central());
      CheckDistance<TreeType, MetricType>(tree, node->Inner());
      CheckDistance<TreeType, MetricType>(tree, node->Outer());
    }
  }
}

template<typename TreeType>
void CheckHRectBound(const TreeType& tree)
{
  if (tree.IsLeaf())
  {
    for (size_t k = 0; k < tree.NumPoints(); k++)
      BOOST_REQUIRE_EQUAL(true,
          tree.Bound().Contains(tree.Dataset().col(tree.Point(k))));
    return;
  }

  for (size_t k = 0; k < tree.NumDescendants(); k++)
  {
    if (!tree.Bound().Contains(tree.Dataset().col(tree.Descendant(k))))
      BOOST_REQUIRE_EQUAL(true,
          tree.Bound().Contains(tree.Dataset().col(tree.Descendant(k))));
  }

  CheckHRectBound(*tree.Central());
  CheckHRectBound(*tree.Inner());
  CheckHRectBound(*tree.Outer());
}

BOOST_AUTO_TEST_CASE(HRectVPTreeDistanceTest)
{
  typedef HRectVPTree<EuclideanDistance, EmptyStatistic, arma::mat> TreeType;
  arma::mat dataset(10, 1000);

  dataset.randn();

  TreeType tree(dataset);
  CheckDistance<TreeType, EuclideanDistance>(tree);
  CheckHRectBound(tree);
}

BOOST_AUTO_TEST_CASE(HRectVPTreeDualTreeTraverserTest)
{
  arma::mat dataset;
  dataset.randu(10, 1000); // 1000 points in 8 dimensions.
  arma::Mat<size_t> neighbors1;
  arma::mat distances1;
  arma::Mat<size_t> neighbors2;
  arma::mat distances2;

  // Nearest neighbor search with the VP tree.
  NeighborSearch<NearestNeighborSort, metric::LMetric<2, true>, arma::mat,
      HRectVPTree> knn1(dataset, false, false);

  knn1.Search(5, neighbors1, distances1);

  // Nearest neighbor search the naive way.
  KNN knn2(dataset, true, true);

  knn2.Search(5, neighbors2, distances2);

  for (size_t i = 0; i < neighbors1.size(); i++)
  {
      BOOST_REQUIRE_EQUAL(neighbors1[i], neighbors2[i]);
      BOOST_REQUIRE_EQUAL(distances1[i], distances2[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END();
