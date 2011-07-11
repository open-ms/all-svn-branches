// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:expandtab
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2011 -- Bastian Blank
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <queue>
#include <boost/shared_ptr.hpp>
#include <boost/unordered/unordered_set.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

namespace OpenMS
{
  /**
   * @brief generic n-dimensional hierarchical clustering
   */
  template <typename PointRef>
  class HierarchicalClustering
  {
    public:
      typedef boost::array<DoubleReal, 2> Point;

     /**
       * @brief Cluster of points.
       *
       * Describes a set of points in the grid.
       */
      typedef typename boost::unordered_multimap<Point, PointRef> Cluster;
      typedef HashGrid<Cluster, 2> Grid;

      Grid grid;

    protected:
      /** @brief Bounding box of cluster. */
      // XXX: wo anders schon definiert
      class BoundingBox
        : public std::pair<Point, Point>
      {
        public:
          BoundingBox(const Point &p)
            : std::pair<Point, Point>(std::make_pair(p, p))
          { }

          Point size() const
          {
            return point_minus(this->second, this->first);
          }

          BoundingBox &operator|=(const BoundingBox &rhs)
          {
            typename Point::iterator lit;
            typename Point::const_iterator rit;

            // Calculate lower bound
            lit = this->first.begin(); rit = rhs.first.begin();
            for (; lit != this->first.end(); ++lit, ++rit) *lit = std::min(*lit, *rit);

            // Calculate upper bound
            lit = this->second.begin(); rit = rhs.second.begin();
            for (; lit != this->second.end(); ++lit, ++rit) *lit = std::max(*lit, *rit);

            return *this;
          }

          BoundingBox operator|(const BoundingBox &rhs) const
          {
            BoundingBox ret(*this);
            ret |= rhs;
            return ret;
          }

          operator Point() const
          {
            // (first + second) / 2
            return point_division(point_plus(this->first, this->second), 2);
          }
      };

      /** @brief Tree node used for clustering. */
      class TreeNode
      {
        public:
          const BoundingBox bbox;
          const boost::shared_ptr<TreeNode> left, right;
          const bool center;

          TreeNode(BoundingBox bbox, bool center)
            : bbox(bbox), center(center)
          { }

          TreeNode(BoundingBox bbox, boost::shared_ptr<TreeNode> &left, boost::shared_ptr<TreeNode> &right)
            : bbox(bbox), left(left), right(right), center(left->center | right->center)
          { }
      };

      /** @brief Distance info used for clustering. */
      // XXX: Use boost::tuple
      // XXX: k-d-tree???
      typedef std::priority_queue<std::pair<DoubleReal, std::pair<boost::shared_ptr<TreeNode>, boost::shared_ptr<TreeNode> > > > DistanceQueue;

    public:
      HierarchicalClustering(const Point &max_delta)
        : grid(max_delta)
      { }

      /**
       * @brief Insert new Point into grid.
       * @param d Point to insert.
       * @param ref Associated caller specified info.
       * @return iterator to inserted cluster.
       */
      typename Grid::local_iterator insertPoint(const Point &d, const PointRef &ref)
      {
        typename Grid::local_iterator it = insertCluster(d);
        it->second.insert(std::make_pair(d, ref));
        return it;
      }

      /**
       * @brief Perform clustering of all existing points.
       */
      void cluster()
      {
        std::cout << "cluster: coord: " << grid.max_key[0] << ":" << grid.max_key[1] << std::endl;
        // Collect coordinates of all active cells
        std::vector<typename Grid::CellIndex> cells;
        for (typename Grid::const_cell_iterator it = grid.cell_begin(); it != grid.cell_end(); ++it)
          cells.push_back(it->first);
        // Cluster each available cell
        for (typename std::vector<typename Grid::CellIndex>::const_iterator it = cells.begin(); it != cells.end(); ++it)
          clusterCell(*it);
      }

    protected:
      /**
       * @brief Insert new Cluster into grid.
       * @param d Point to insert.
       * @return iterator to inserted cluster.
       */
      typename Grid::local_iterator insertCluster(const Point &d)
      {
        return grid.insert(std::make_pair(d, Cluster()));
      }

      typedef std::map<typename Grid::CellIndex, std::pair<const typename Grid::Cell *, bool> > ClusterCells;
      typedef boost::shared_ptr<TreeNode> ClusterTree;

      /**
       * @brief Perform clustering for each available cell.
       * @param p Cell coordinate to cluster.
       */
      void clusterCell(const typename Grid::CellIndex &p);

      void clusterCellCollect(UInt dim, typename Grid::CellIndex cur, ClusterCells cells, bool center);

      void clusterCellReaddCluster(const ClusterTree &tree, Cluster &cluster)
      {
        if (tree->left && tree->right)
        {
          clusterCellReaddCluster(tree->left, cluster);
          clusterCellReaddCluster(tree->right, cluster);
        }
        else
        {
          // XXX: PointRef
          cluster.insert(std::make_pair(tree->bbox.first, PointRef()));
        }
      }

      void clusterCellReaddPoint(const ClusterTree &tree)
      {
        if (tree->left && tree->right)
        {
          clusterCellReaddPoint(tree->left);
          clusterCellReaddPoint(tree->right);
        }
        else
        {
          // XXX: PointRef
          insertPoint(tree->bbox.first, PointRef());
        }
      }

      // XXX: Convert to operator
      static Point point_plus(const Point &lhs, const Point &rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit + *rit;
        return ret;
      }

      // XXX: Convert to operator
      static Point point_minus(const Point &lhs, const Point &rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit - *rit;
        return ret;
      }

      // XXX: Convert to operator
      static Point point_multiplication(const Point &lhs, const Point &rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit * *rit;
        return ret;
      }

      // XXX: Convert to operator
      static Point point_division(const Point &lhs, const DoubleReal &rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin();
        for (; it != ret.end(); ++it, ++lit) *it = *lit / rhs;
        return ret;
      }

      // XXX: Convert to operator
      static Point point_division(const Point &lhs, const Point &rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; it != ret.end(); ++it, ++lit, ++rit) *it = *lit / *rit;
        return ret;
      }

      // XXX: Convert to operator
      static bool point_greater(const Point &lhs, const Point &rhs)
      {
        UInt ret = 0;
        typename Point::const_iterator lit = lhs.begin(), rit = rhs.begin();
        for (; lit != lhs.end(); ++lit, ++rit) ret += *lit > *rit;
        return ret;
      }

      static DoubleReal point_distance(const Point &lhs, const Point &rhs)
      {
        DoubleReal ret = 0;
        Point p = point_elem_power(point_minus(lhs, rhs), 2);
        typename Point::const_iterator it = p.begin();
        for (; it != p.end(); ++it) ret += *it;
        return std::sqrt(ret);
      }

      static Point point_elem_power(const Point &lhs, const UInt rhs)
      {
        Point ret;
        typename Point::iterator it = ret.begin();
        typename Point::const_iterator lit = lhs.begin();
        for (; it != ret.end(); ++it, ++lit) *it = std::pow(*lit, rhs);
        return ret;
      }

  };

  template <typename I>
  void HierarchicalClustering<I>::clusterCell(const typename Grid::CellIndex &cur)
  {
    typedef boost::unordered_set<ClusterTree> LocalTrees;

    ClusterCells cells;
    LocalTrees trees;
    DistanceQueue dists;

    // Collect all cells we need
    std::cout << "ping: coord: " << cur[0] << ":" << cur[1] << std::endl;
    try { cells.insert(std::make_pair(cur, std::make_pair(&grid.cell_at(cur), true))); }
    catch (std::out_of_range &) { return; }
    clusterCellCollect(1, cur, cells, true);

    // Collect existing points per cell
    std::cout << "ping" << std::endl;
    for (typename ClusterCells::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it)
    {
      const typename Grid::Cell &cell_cur = *cell_it->second.first;
      const bool &cell_center = cell_it->second.second;

      // Per cluster
      for (typename Grid::const_local_iterator cluster_it = cell_cur.begin(); cluster_it != cell_cur.end(); ++cluster_it)
      {
        // Per point
        for (typename Cluster::const_iterator point_it = cluster_it->second.begin(); point_it != cluster_it->second.end(); ++point_it)
        {
          boost::shared_ptr<TreeNode> tree(new TreeNode(point_it->first, cell_center));

          // Generate distance to every existing tree
          Point new_normpoint = point_multiplication(tree->bbox, grid.max_delta);
          for (typename LocalTrees::const_iterator it = trees.begin(); it != trees.end(); ++it)
          {
            Point old_normpoint = point_multiplication((*it)->bbox, grid.max_delta);
            DoubleReal dist = point_distance(new_normpoint, old_normpoint);
            dists.push(std::make_pair(dist, std::make_pair(tree, *it)));
          }

          trees.insert(tree);
        }
      }
    }

    // Join points
    std::cout << "ping: size: " << trees.size() << ", " << dists.size() << std::endl;
    while (!dists.empty())
    {
      const typename DistanceQueue::value_type &cur_dist = dists.top();
      boost::shared_ptr<TreeNode> tree_left(cur_dist.second.first), tree_right(cur_dist.second.second);
      dists.pop();
 
      // Chck if both trees are still available
      if (trees.count(tree_left) + trees.count(tree_right) < 2)
      {
        continue;
      }

      // Calculate new bounding box.
      BoundingBox bbox = tree_left->bbox | tree_right->bbox;
      boost::shared_ptr<TreeNode> tree(new TreeNode(bbox, tree_left, tree_right));

      trees.erase(tree_left);
      trees.erase(tree_right);

      if (!point_greater(bbox.size(), grid.max_delta))
      {
        // Generate distance to every existing tree
        // XXX: De-duplicate
        Point new_normpoint = point_multiplication(tree->bbox, grid.max_delta);
        for (typename LocalTrees::const_iterator it = trees.begin(); it != trees.end(); ++it)
        {
          Point old_normpoint = point_multiplication((*it)->bbox, grid.max_delta);
          DoubleReal dist = point_distance(new_normpoint, old_normpoint);
          dists.push(std::make_pair(dist, std::make_pair(tree, *it)));
        }
      }

      trees.insert(tree);
    }

    // Clear cells
    std::cout << "ping" << std::endl;
    for (typename ClusterCells::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it)
    {
      grid.cell_clear(cell_it->first);
    }

    // Add current data to grid
    std::cout << "ping: size: " << trees.size() << ", " << dists.size() << std::endl;
    for (typename LocalTrees::iterator tree_it = trees.begin(); tree_it != trees.end(); ++tree_it)
    {
      // We got a finished tree with a point in the center, add cluster
      if ((**tree_it).center)
      {
        Cluster &cluster = insertCluster((**tree_it).bbox)->second;
        clusterCellReaddCluster(*tree_it, cluster);
      }
      // We got a finished tree but no point in the center, readd as single points
      else
      {
        clusterCellReaddPoint(*tree_it);
      }
    }
  }

  // XXX: no recursion
  template <typename I>
  void HierarchicalClustering<I>::clusterCellCollect(UInt dim, typename Grid::CellIndex cur, ClusterCells cells, bool center)
  {
    if (dim)
    {
      clusterCellCollect(dim - 1, cur, cells, center);
      cur[dim] -= 1;
      clusterCellCollect(dim - 1, cur, cells, false);
      cur[dim] += 2;
      clusterCellCollect(dim - 1, cur, cells, false);
    }
    else if (!center)
    {
      try
      {
        cells.insert(std::make_pair(cur, std::make_pair(&grid.cell_at(cur), false)));
      }
      catch (std::out_of_range &)
      { }
    }
  }
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
