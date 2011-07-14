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
#include <boost/unordered/unordered_set.hpp>

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/CONCEPT/Types.h>

#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

namespace OpenMS
{
  /**
   * @brief Generic 2-dimensional hierarchical clustering.
   *
   * @tparam PointRef Caller specified referenced associated with every point.
   * Needs to be default constructible.
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

          /** @brief Intersection of bounding box. */
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

          /** @brief Intersection of bounding box. */
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
          const Point normcoord;
          TreeNode *left, *right;
          const bool center;
          const PointRef ref;

          TreeNode(const BoundingBox &bbox, const Point &normcoord, const PointRef &ref, bool center)
            : bbox(bbox), normcoord(normcoord), left(0), right(0), center(center), ref(ref)
          { }

          TreeNode(const BoundingBox &bbox, const Point &normcoord, TreeNode *left, TreeNode *right)
            : bbox(bbox), normcoord(normcoord), left(left), right(right), center(left->center && right->center), ref(PointRef())
          { }
      };

      /** @brief Distance info used for clustering. */
      class DistanceInfo
      {
        public:
          DoubleReal distance;
          TreeNode *left, *right;

          DistanceInfo(const DoubleReal &distance, TreeNode *left, TreeNode *right)
            : distance(distance), left(left), right(right)
          { }

          bool operator<(const DistanceInfo &rhs) const
          {
            return (distance < rhs.distance ||
                left < rhs.left ||
                right < rhs.right);
          }
      };

      /** @brief Distance queue used for clustering. */
      // XXX: k-d-tree???
      typedef std::priority_queue<DistanceInfo> DistanceQueue;

    public:
      /**
       * @brief Constructor
       * @param max_delta Max size of cluster
       */
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

      typedef std::map<typename Grid::CellIndex, std::pair<typename Grid::Cell *, bool> > ClusterCells;

      /**
       * @brief Perform clustering on given cell.
       * @param p Cell index.
       */
      void clusterCell(const typename Grid::CellIndex &p);

      /**
       * @brief Collect cells used for given cell.
       * This function collects all cells in a 5x5 array.
       * @param cur Cell index.
       * @param cells List of cells to be used.
       */
      void clusterCellCollect(typename Grid::CellIndex cur, ClusterCells &cells);

      /**
       * @brief Collect one cell.
       * @param cur Cell index.
       * @param cells List of cells.
       * @param center Is the given cell in the center.
       * @oaram ignore_missing Defines if non-existant errors should be ignored.
       */
      void clusterCellCollectOne(const typename Grid::CellIndex &cur, ClusterCells &cells, bool center = false, bool ignore_missing = true)
      {
        try
        {
          cells.insert(std::make_pair(cur, std::make_pair(&grid.cell_at(cur), center)));
        }
        catch (std::out_of_range &)
        {
          if (!ignore_missing) throw;
        }
      }

      /**
       * @brief Recursivly readd the points of a finished cluster.
       * All points are saved in the leafs of the tree.
       * @param tree The tree
       * @param cluster The cluster
       */
      void clusterCellReaddCluster(const TreeNode *tree, Cluster &cluster)
      {
        if (tree->left && tree->right)
        {
          clusterCellReaddCluster(tree->left, cluster);
          clusterCellReaddCluster(tree->right, cluster);
          delete tree->left;
          delete tree->right;
        }
        else
        {
          cluster.insert(std::make_pair(tree->bbox.first, tree->ref));
        }
      }

      /**
       * @brief Recursively readd the points of an unfinished cluster back to the grid.
       * All points are saved in the leafs of the tree.
       * @param tree The tree
       */
      void clusterCellReaddPoint(const TreeNode *tree)
      {
        if (tree->left && tree->right)
        {
          clusterCellReaddPoint(tree->left);
          clusterCellReaddPoint(tree->right);
        }
        else
        {
          insertPoint(tree->bbox.first, tree->ref);
        }
        delete tree->left;
        delete tree->right;
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
        Point p = point_minus(lhs, rhs);
        typename Point::const_iterator it = p.begin();
        for (; it != p.end(); ++it) ret += std::pow(*it, 2.);
        return std::sqrt(ret);
      }
  };

  template <typename I>
  void HierarchicalClustering<I>::clusterCell(const typename Grid::CellIndex &cur)
  {
    typedef boost::unordered_set<TreeNode *> LocalTrees;

    ClusterCells cells;
    LocalTrees trees;
    DistanceQueue dists;

    // Collect all cells we need
    std::cout << "start: coord: " << cur[0] << ":" << cur[1] << '\n';
    try
    {
      clusterCellCollect(cur, cells);
    }
    catch (std::out_of_range &)
    { return; }

    // Collect and remove existing points from cells
    std::cout << "number of cells: " << cells.size() << '\n';
    int rounds = 0;
    for (typename ClusterCells::iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it)
    {
      typename Grid::Cell &cell_cur = *cell_it->second.first;
      const bool &cell_center = cell_it->second.second;

      // Iterate per cluster
      typename Grid::local_iterator cluster_tmp_it = cell_cur.begin();
      while (cluster_tmp_it != cell_cur.end())
      {
        typename Grid::local_iterator cluster_it = cluster_tmp_it;
        ++cluster_tmp_it;

        // Not yet a cluster
        if (cluster_it->second.size() == 1)
        {
          // Per point
          for (typename Cluster::const_iterator point_it = cluster_it->second.begin(); point_it != cluster_it->second.end(); ++point_it)
          {
            rounds++;
            if (!(rounds % 1000)) std::cout << "  ping:" << rounds << '\n';
            const BoundingBox bbox = point_it->first;
            const Point normcoord = point_division(point_it->first, grid.max_delta);
            TreeNode *tree(new TreeNode(bbox, normcoord, point_it->second, cell_center));

            // Generate distance to every existing tree
            for (typename LocalTrees::const_iterator it = trees.begin(); it != trees.end(); ++it)
            {
              DoubleReal dist = point_distance(tree->normcoord, (*it)->normcoord);
              dists.push(DistanceInfo(dist, tree, *it));
            }

            trees.insert(tree);
          }

          cell_cur.erase(cluster_it);
        }
      }
    }

    // Join points
    std::cout << "initial trees: size: " << trees.size() << ", " << dists.size() << '\n';
    rounds = 0;
    while (!dists.empty())
    {
      rounds++;
      if (!(rounds % 1000000)) std::cout << "  ping:" << rounds << '\n';
      const typename DistanceQueue::value_type &cur_dist = dists.top();
      TreeNode *tree_left(cur_dist.left), *tree_right(cur_dist.right);
      dists.pop();
 
      // Chck if both trees are still available
      if (trees.count(tree_left) + trees.count(tree_right) < 2)
      {
        continue;
      }

      const BoundingBox bbox = tree_left->bbox | tree_right->bbox;

      if (!point_greater(bbox.size(), grid.max_delta))
      {
        const Point normcoord = point_division(bbox, grid.max_delta);
        TreeNode *tree(new TreeNode(bbox, normcoord, tree_left, tree_right));
        trees.erase(tree_left);
        trees.erase(tree_right);

        // Generate distance to every existing tree
        // XXX: De-duplicate
        for (typename LocalTrees::const_iterator it = trees.begin(); it != trees.end(); ++it)
        {
          DoubleReal dist = point_distance(tree->normcoord, (*it)->normcoord);
          dists.push(DistanceInfo(dist, tree, *it));
        }

        trees.insert(tree);
      }
    }

    // Add current data to grid
    std::cout << "late trees: size: " << trees.size() << ", " << dists.size() << "; rounds: " << rounds << '\n';
    for (typename LocalTrees::iterator tree_it = trees.begin(); tree_it != trees.end(); ++tree_it)
    {
      // We got a finished tree with all points in the center, add cluster
      if ((**tree_it).center)
      {
        Cluster &cluster = insertCluster((**tree_it).bbox)->second;
        clusterCellReaddCluster(*tree_it, cluster);
      }
      // We got a finished tree but not all points in the center, readd as single points
      else
      {
        clusterCellReaddPoint(*tree_it);
      }
      delete *tree_it;
    }

    std::cout << "end\n";
  }

  // XXX: check if 2x2 center and 4x4 is sufficient
  template <typename I>
  void HierarchicalClustering<I>::clusterCellCollect(typename Grid::CellIndex base, ClusterCells &cells)
  {
    // (0, 0)
    clusterCellCollectOne(base, cells, true, false);

    typename Grid::CellIndex cur = base;
    cur[0] -= 1;
    // (-1, -1)
    cur[1] -= 1; clusterCellCollectOne(cur, cells);
    // (-1, 0)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (-1, 1)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (-1, 2)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (-1, 3)
    cur[1] += 1; clusterCellCollectOne(cur, cells);

    cur = base;
    // (0, -1)
    cur[1] -= 1; clusterCellCollectOne(cur, cells);
    // (0, 0)
    cur[1] += 1;
    // (0, 1)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (0, 2)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (0, 3)
    cur[1] += 1; clusterCellCollectOne(cur, cells);

    cur = base; cur[0] += 1;
    // (1, -1)
    cur[1] -= 1; clusterCellCollectOne(cur, cells);
    // (1, 0)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (1, 1)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (1, 2)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (1, 3)
    cur[1] += 1; clusterCellCollectOne(cur, cells);

    cur = base; cur[0] += 2;
    // (2, -1)
    cur[1] -= 1; clusterCellCollectOne(cur, cells);
    // (2, 0)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (2, 1)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (2, 2)
    cur[1] += 1; clusterCellCollectOne(cur, cells, true);
    // (2, 3)
    cur[1] += 1; clusterCellCollectOne(cur, cells);

    cur = base; cur[0] += 3;
    // (3, -1)
    cur[1] -= 1; clusterCellCollectOne(cur, cells);
    // (3, 0)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (3, 1)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (3, 2)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
    // (3, 3)
    cur[1] += 1; clusterCellCollectOne(cur, cells);
  }
}

#endif /* OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H */
