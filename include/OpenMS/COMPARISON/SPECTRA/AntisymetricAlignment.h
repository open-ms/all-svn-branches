// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_ANTISYMETRICALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_ANTISYMETRICALIGNMENT_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#define ALIGNMENT_DEBUG
#undef  ALIGNMENT_DEBUG

namespace OpenMS
{

	//~ warning: unused parameter 'parent_mass2'
	//~ warning: unused parameter 'common_scores'
	//~ warning: unused parameter 'common2_scores'

	struct AntisymetricDP
	{

				enum celltype {cell_prefix1, cell_suffix1, cell_prefix2, cell_suffix2,
							cell_prefix1_L, cell_suffix1_R, cell_prefix2_L, cell_suffix2_R,
							cell_D1, cell_D2, cell_D3,
							cell_M1_L, cell_M1_R, cell_M2_L, cell_M2_R,
							cell_M3_L, cell_M3_R, INVALID};

		std::pair<double, std::pair< std::vector<int>,std::vector<int> > > calculateAlignementPath(double parent_mass1,  double parent_mass2, std::vector<double> peaks, std::vector<double> peaks2, std::vector<int> & common, std::vector<int> & common2, std::vector<double> & common_scores, std::vector<double> & common2_scores, std::vector<int> & prev, std::vector<int> & next, std::vector<int> & prev2, std::vector<int> & next2, std::vector<std::vector<int> > & left_jumps, std::vector<std::vector<int> > & right_jumps, std::vector<std::vector<int> > & left_jumps2, std::vector<std::vector<int> > &right_jumps2, std::vector<std::vector<int> > & left_neighbors, std::vector<std::vector<int> > & right_neighbors, Real peak_tolerance, double same_vertex_penalty, double ptm_penalty)
		{

			const DoubleReal infinity = std::numeric_limits<DoubleReal>::max();

			int n0 = common.size();
			int N = n0-1;
			int n = (n0+1)/2;

			std::vector<double> prefix1(n, -infinity);
			std::vector<std::vector<double> > suffix1(n);
			std::vector<std::vector<double> > prefix2(n);
			std::vector<double> suffix2(n, -infinity);

			std::vector<double> prefix1_L(n, -infinity);
			std::vector<double> suffix1max(n, -infinity);
			std::vector<double> suffix1_R(n, -infinity);
			std::vector<double> prefix2max(n, -infinity);
			std::vector<double> prefix2_L(n, -infinity);
			std::vector<double> suffix2_R(n, -infinity);

			std::vector<std::vector<double> > D1(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<std::vector<double> > > D2(n,std::vector<std::vector<double> >(n));
			std::vector<std::vector<std::vector<double> > > D3(n,std::vector<std::vector<double> >(n));
			std::vector<std::vector<double> > D2max(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<double> > D3max(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<double> > M1_L(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<double> > M1_R(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<double> > M2_L(n,std::vector<double>(n, -infinity));
			std::vector<std::vector<std::vector<double> > > M2_R(n,std::vector<std::vector<double> >(n));
			std::vector<std::vector<std::vector<double> > > M3_L(n,std::vector<std::vector<double> >(n));
			std::vector<std::vector<double> > M3_R(n,std::vector<double>(n, -infinity));

			// compute prefix1
			for(int i = 0; i < n; ++i)
			{
				if(common[i] != -1)
				{
					double prev_score = 0.0;
					for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
					{
						prev_score = std::max(prev_score, prefix1[*it]);
					}
					int i2 = prev[i];
					if (i2 != -1)
					{
						prev_score = std::max(prev_score, prefix1_L[i2]);
					}
					prefix1[i] = common_scores[i] + prev_score;
				}

				if(i != 0)
				{
					prefix1_L[i] = std::max(prefix1_L[i-1], prefix1[i]);
				}
				else
				{
					prefix1_L[i] = prefix1[i];
				}
			}

			// compute suffix1
			for(int j = 0; j < n; ++j)
			{
				int l = right_neighbors[N-j].size();
				suffix1[j].resize(l+1, -infinity);

				if (common[N-j] != -1)
				{ // s = 0
					double next_score = 0.0;
					for(std::vector<int>::const_iterator it = right_jumps[N-j].begin();
					it != right_jumps[N-j].end(); ++it)
					{
						next_score = std::max(next_score, suffix2[N-*it]+ptm_penalty);
					}
					int j2 = next[N-j];
					if(j2 != -1)
					{
						next_score = std::max(next_score, suffix2_R[N-j2]+ptm_penalty);
					}
					for(std::vector<int>::const_iterator it = right_jumps2[N-j].begin(); it != right_jumps2[N-j].end(); ++it)
					{
						next_score = std::max(next_score, suffix1max[N-*it]);
					}
					j2 = next2[N-j];
					if(j2 != -1)
					{
						next_score = std::max(next_score, suffix1_R[N-j2]);
					}
					suffix1[j][0] = common_scores[N-j] + next_score;

					// s > 0
					for(int s = 0; s < l; ++s)
					{
						int j2 = right_neighbors[N-j][s];
						if (common[j2] == -1)
						{
							continue;
						}
						double penalty = 0.0;
						double y = peaks2[common[N-j]]+peaks2[common[j2]];
						if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
						penalty = same_vertex_penalty;
						double next_score = suffix1max[N-j2];
						suffix1[j][s+1] = common_scores[N-j] + next_score + penalty;
					}
				}

				suffix1max[j] = *std::max_element(suffix1[j].begin(), suffix1[j].end());
				if(j != 0)
				{
					suffix1_R[j] = std::max(suffix1_R[j-1], suffix1max[j]);
				}
				else
				{
					suffix1_R[j] = suffix1max[j];
				}
			}

			// compute prefix2
			for(int i = 0; i < n; ++i)
			{
				int l = left_neighbors[i].size();
				prefix2[i].resize(l+1 ,-infinity);

				if (common2[i] != -1)
				{// s = 0
					double prev_score = 0.0;
					for (std::vector<int>::const_iterator it = left_jumps[i].begin();it != left_jumps[i].end(); ++it)
					{
						prev_score = std::max(prev_score, prefix1[*it]+ptm_penalty);
					}
					int i2 = prev[i];
					if (i2 != -1)
					{
						prev_score = std::max(prev_score, prefix1_L[i2]+ptm_penalty);
					}
					for(std::vector<int>::const_iterator it = left_jumps2[i].begin(); it != left_jumps2[i].end(); ++it)
					{
						prev_score = std::max(prev_score, prefix2max[*it]);
					}
					i2 = prev2[i];
					if (i2 != -1)
					{
						prev_score = std::max(prev_score, prefix2_L[i2]);
					}
					prefix2[i][0] = common2_scores[i] + prev_score;

					// s > 0
					for (int s = 0; s < l; ++s)
					{
						int i2 = left_neighbors[i][s];
						if(common2[i2] == -1)
						{
							continue;
						}
						double penalty = 0.0;
						double y = peaks2[common2[i]]+peaks2[common2[i2]];
						if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
						{
							penalty = same_vertex_penalty;
						}
						prev_score = prefix2max[i2];
						prefix2[i][s+1] = common2_scores[i] + prev_score + penalty;
					}
				}

				prefix2max[i] = *std::max_element(prefix2[i].begin(), prefix2[i].end());
				if (i != 0)
				{
					prefix2_L[i] = std::max(prefix2_L[i-1], prefix2max[i]);
				}
				else
				{
					prefix2_L[i] = prefix2max[i];
				}
			}

			// compute suffix2
			for(int j = 0; j < n; ++j)
			{
				if(common2[N-j] != -1)
				{
					double next_score = 0.0;
					for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
					{
						next_score = std::max(next_score, suffix2[N-*it]);
					}
					int j2 = next[N-j];
					if(j2 != -1)
					{
						next_score = std::max(next_score, suffix2_R[N-j2]);
					}
					suffix2[j] = common2_scores[N-j] + next_score;
				}
				if(j != 0)
				{
					suffix2_R[j] = std::max(suffix2_R[j-1], suffix2[j]);
				}
				else
				{
					suffix2_R[j] = suffix2[j];
				}
			}

			// compute D1/M1_R
			for (int i = 0; i < n; ++i)
			{
				if (common[i] == -1)
				{
					if (i != 0)
					{
						M1_L[i] = M1_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(common2[N-j] == -1)
					{
						if(i != 0)
						{
							M1_L[i][j] = M1_L[i-1][j];
						}
						if(j != 0)
						{
							M1_R[i][j] = M1_R[i][j-1];
						}
						continue;
					}

					double penalty = 0.0;
					double x = peaks[i]+peaks[N-j];
					double y = peaks2[common[i]]+peaks2[common2[N-j]];
					if(fabs(x-(parent_mass1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = same_vertex_penalty;
					}

					if(i >= j)
					{
						double prev_score = suffix2[j]+ptm_penalty;
						// If we choose suffix2[j], then we pay the PTM penalty - the penalty is paid only once!
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, D1[*it][j]);
						}
						int i2 = prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, M1_L[i2][j]);
						}
						D1[i][j] = common_scores[i] + prev_score + penalty;
					}
					else
					{
						double next_score = prefix1[i]+ptm_penalty;
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, D1[i][N-*it]);
						}
						int j2 = next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, M1_R[i][N-j2]);
						}
						D1[i][j] = common2_scores[N-j] + next_score + penalty;
					}

					// Compute M1_L
					if(i != 0)
					{
						M1_L[i][j] = std::max(M1_L[i-1][j], D1[i][j]);
					}
					else
					{
						M1_L[i][j] = D1[i][j];
					}

					// Compute M1_R
					if(j != 0)
					{
						M1_R[i][j] = std::max(M1_R[i][j-1], D1[i][j]);
					}
					else
					{
						M1_R[i][j] = D1[i][j];
					}
				}
			}

			// compute D2/D2max/M2_L
			for (int i = 0; i < n; ++i)
			{
				int l = left_neighbors[i].size();
				for(int j = 0; j < n; ++j)
				{
					D2[i][j].resize(l+1, -infinity);
					M2_R[i][j].resize(l+1, -infinity);
				}
				if(common2[i] == -1)
				{
					if (i > 0)
					{
						M2_L[i] = M2_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(common2[N-j] == -1)
					{
						if(i != 0)
						{
							M2_L[i][j] = M2_L[i-1][j];
						}
						if(j != 0)
						{
							M2_R[i][j] = M2_R[i][j-1];
						}
						continue;
					}

					double penalty = 0.0;
					double x = peaks[i]+peaks[N-j];
					double y = peaks2[common2[i]]+peaks2[common2[N-j]];
					if(fabs(x-(parent_mass1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = same_vertex_penalty;
					}

					if (i > j)
					{// s = 0
						double prev_score = suffix2[j];
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, D1[*it][j]);
						}
						int i2 = prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, M1_L[i2][j]);
						}
						for(std::vector<int>::const_iterator it = left_jumps2[i].begin(); it != left_jumps2[i].end(); ++it)
						{
							prev_score = std::max(prev_score, D2max[*it][j]);
						}
						i2 = prev2[i];
						if (i2 != -1)
						{
							prev_score = std::max(prev_score, M2_L[i2][j]);
						}
						D2[i][j][0] = common2_scores[i] + prev_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int i2 = left_neighbors[i][s];
							if(common2[i2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = peaks2[common2[i]]+peaks2[common2[i2]];
							if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += same_vertex_penalty;
							}
							prev_score = D2max[i2][j];
							D2[i][j][s+1] = common2_scores[i] + prev_score + penalty2;
						}
					}
					else
					{// i <= j
						// s = 0
						double next_score = prefix2[i][0];
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, D2[i][N-*it][0]);
						}
						int j2 = next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, M2_R[i][N-j2][0]);
						}
						D2[i][j][0] = common2_scores[N-j] + next_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int i2 = left_neighbors[i][s];
							if(common2[i2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = peaks2[common2[i2]]+peaks2[common2[N-j]];
							if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += same_vertex_penalty;
							}
							double next_score = prefix2[i][s+1];
							for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
							{
								next_score = std::max(next_score, D2[i][N-*it][s+1]);
							}
							int j2 = next[N-j];
							if(j2 != -1)
							{
								next_score = std::max(next_score, M2_R[i][N-j2][s+1]);
							}
							D2[i][j][s+1] = common2_scores[N-j] + next_score + penalty2;
						}
					}

					// Compute D2max
					D2max[i][j] = *std::max_element(D2[i][j].begin(), D2[i][j].end());

					// Compute M2_L
					if(i != 0)
					{
						M2_L[i][j] = std::max(M2_L[i-1][j], D2max[i][j]);
					}
					else
					{
						M2_L[i][j] = D2max[i][j];
					}

					// compute M2_R
					if(j != 0)
					{
						for(int s = 0; s < l+1; ++s)
						{
							M2_R[i][j][s] = std::max(M2_R[i][j-1][s], D2[i][j][s]);
						}
					}
					else
					{
						M2_R[i][j] = D2[i][j];
					}
				}
			}

			// compute D3/M3_R
			for(int i = 0; i < n; ++i)
			{
				for(int j = 0; j < n; ++j)
				{
					int l = right_neighbors[N-j].size();
					D3[i][j].resize(l+1, -infinity);
					M3_L[i][j].resize(l+1, -infinity);
				}
				if(common[i] == -1)
				{
					if(i > 0)
					{
						M3_L[i] = M3_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(common[N-j] == -1)
					{
						if(i != 0)
						{
							M3_L[i][j] = M3_L[i-1][j];
						}
						if(j != 0)
						{
							M3_R[i][j] = M3_R[i][j-1];
						}
						continue;
					}

					int l = right_neighbors[N-j].size();
					double penalty = 0.0;
					double x = peaks[i]+peaks[N-j];
					double y = peaks2[common[i]]+peaks2[common[N-j]];
					if(fabs(x-(parent_mass1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = same_vertex_penalty;
					}

					if(i >= j)
					{ // s = 0
						double prev_score = suffix1[j][0];
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, D3[*it][j][0]);
						}
						int i2 = prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, M3_L[i2][j][0]);
						}
						D3[i][j][0] = common_scores[i] + prev_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int j2 = right_neighbors[N-j][s];
							if(common[j2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = peaks2[common[i]]+peaks2[common[j2]];
							if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += same_vertex_penalty;
							}
							double prev_score = suffix1[j][s+1];
							for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
							{
								prev_score = std::max(prev_score, D3[*it][j][s+1]);
							}
							int i2 = prev[i];
							if(i2 != -1)
							{
								prev_score = std::max(prev_score, M3_L[i2][j][s+1]);
							}
							D3[i][j][s+1] = common_scores[i] + prev_score + penalty2;
						}
					}
					else
					{ // s = 0
						double next_score = prefix1[i];
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, D1[i][N-*it]);
						}
						int j2 = next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, M1_R[i][N-j2]);
						}
						for(std::vector<int>::const_iterator it = right_jumps2[N-j].begin(); it != right_jumps2[N-j].end(); ++it)
						{
							next_score = std::max(next_score, D3max[i][N-*it]);
						}
						j2 = next2[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, M3_R[i][N-j2]);
						}
						D3[i][j][0] = common_scores[N-j] + next_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int j2 = right_neighbors[N-j][s];
							if(common[j2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = peaks2[common[N-j]]+peaks2[common[j2]];
							if(fabs(y-(parent_mass2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += same_vertex_penalty;
							}
							double next_score = D3max[i][N-j2];
							D3[i][j][s+1] = common_scores[N-j] + next_score + penalty2;
						}
					}

					// Compute D3max
					D3max[i][j] = *std::max_element(D3[i][j].begin(), D3[i][j].end());

					// Compute M3_L
					if (i != 0)
					{
						for(int s = 0; s < l+1; ++s)
						{
							M3_L[i][j][s] = std::max(M3_L[i-1][j][s], D3[i][j][s]);
						}
					}
					else
					{
						M3_L[i][j] = D3[i][j];
					}

					// Compute M3_R
					if(j != 0)
					{
						M3_R[i][j] = std::max(M3_R[i][j-1], D3max[i][j]);
					}
					else
					{
						M3_R[i][j] = D3max[i][j];
					}
				}
			}

			// Find best score
			double best_score = 0.0;
			int best_i = -1;
			int best_j = -1;
			int best_s = -1;
			celltype best_t = INVALID;
			for(int i = 0; i < n; ++i)
			{
				for(std::vector<int>::const_iterator it = right_jumps[i].begin(); it != right_jumps[i].end(); ++it)
				{
					int j = N-*it;
					if(j <= n-1)
					{
						if(best_score < D1[i][j])
						{
							best_score = D1[i][j];
							best_i = i;
							best_j = j;
							best_s = 0;
							best_t = cell_D1;
						}
						if(best_score < D2max[i][j])
						{
							best_score = D2max[i][j];
							best_i = i;
							best_j = j;
							best_s = (std::max_element(D2[i][j].begin(),D2[i][j].end())-D2[i][j].begin());
							best_t = cell_D2;
						}
						if(best_score < D3max[i][j])
						{
							best_score = D3max[i][j];
							best_i = i;
							best_j = j;
							best_s = (std::max_element(D3[i][j].begin(),D3[i][j].end())-D3[i][j].begin());
							best_t = cell_D3;
						}
					}
				}

				int j0 = next[i];
				if(j0 != -1)
				{
					int j = std::min(n-1, N-j0);
					if(best_score < M1_R[i][j])
					{
						best_score = M1_R[i][j];
						best_i = i;
						best_j = j;
						best_s = 0;
						best_t = cell_M1_R;
					}

					double tmp = *std::max_element(M2_R[i][j].begin(), M2_R[i][j].end());
					if (best_score < tmp)
					{
						best_score = tmp;
						best_i = i;
						best_j = j;
						best_s = (std::max_element(M2_R[i][j].begin(),M2_R[i][j].end())-M2_R[i][j].begin());
						best_t = cell_M2_R;
					}

					if (best_score < M3_R[i][j])
					{
						best_score = M3_R[i][j];
						best_i = i;
						best_j = j;
						best_s = 0;
						best_t = cell_M3_R;
					}
				}

				if(best_score < prefix1[i])
				{
					best_score = prefix1[i];
					best_i = i;
					best_j = 0;
					best_s = 0;
					best_t = cell_prefix1;
				}

				if(best_score < suffix1max[i])
				{
					best_score = suffix1max[i];
					best_j = i;
					best_i = 0;
					best_s = (std::max_element(suffix1[i].begin(),suffix1[i].end())-suffix1[i].begin());
					best_t = cell_suffix1;
				}

				if(best_score < prefix2max[i])
				{
					best_score = prefix2max[i];
					best_i = i;
					best_j = 0;
					best_s = (std::max_element(prefix2[i].begin(),prefix2[i].end())-prefix2[i].begin());
					best_t = cell_prefix2;
				}

				if(best_score < suffix2[i])
				{
					best_score = suffix2[i];
					best_j = i;
					best_i = 0;
					best_s = 0;
					best_t = cell_suffix2;
				}
			}

			return traceback(parent_mass1, ptm_penalty, prev, next, prev2, next2, left_jumps, right_jumps, left_jumps2, right_jumps2, left_neighbors, right_neighbors, best_i, best_j, best_s, best_t, best_score, D1, D2, D3, D2max, D3max, M1_L, M1_R, M2_L, M2_R, M3_L, M3_R, prefix1, suffix1, prefix2, suffix2, prefix1_L, suffix1max, suffix1_R, prefix2max, prefix2_L, suffix2_R);
		}
		///

		std::pair<double, std::pair< std::vector<int>,std::vector<int> > > traceback(double parent_mass1, double ptm_penalty, std::vector<int> & prev, std::vector<int> & next, std::vector<int> & prev2, std::vector<int> & next2, std::vector<std::vector<int> > & left_jumps, std::vector<std::vector<int> > & right_jumps, std::vector<std::vector<int> > & left_jumps2, std::vector<std::vector<int> > & right_jumps2, std::vector<std::vector<int> > & left_neighbors, std::vector<std::vector<int> > & right_neighbors, int best_i, int best_j, int best_s, celltype best_t, double best_score, std::vector<std::vector<double> > & D1, std::vector<std::vector<std::vector<double> > > & D2, std::vector<std::vector<std::vector<double> > > & D3, std::vector<std::vector<double> > & D2max, std::vector<std::vector<double> > & D3max, std::vector<std::vector<double> > & M1_L, std::vector<std::vector<double> > & M1_R, std::vector<std::vector<double> > & M2_L, std::vector<std::vector<std::vector<double> > > & M2_R, std::vector<std::vector<std::vector<double> > > & M3_L, std::vector<std::vector<double> > & M3_R, std::vector<double> & prefix1, std::vector<std::vector<double> > & suffix1, std::vector<std::vector<double> > & prefix2, std::vector<double> & suffix2, std::vector<double> & prefix1_L, std::vector<double> & suffix1max, std::vector<double> & suffix1_R, std::vector<double> & prefix2max, std::vector<double> & prefix2_L, std::vector<double> & suffix2_R)
		{
			int n0 = prev.size();
			int N = n0-1;

			std::vector<int> path1;
			std::vector<int> path2;

			if (best_i == -1)
			{
				return std::pair<double, std::pair< std::vector<int>,std::vector<int> > >(best_score/parent_mass1, std::pair<std::vector<int>,std::vector<int> >(path1, path2));
			}

			int i = best_i;
			int j = best_j;
			int s = best_s;
			celltype t = best_t;
			while(1)
			{
				if(t == cell_prefix1)
				{
					path1.insert(path1.begin(), i);
					double prev_score = 0.0;
					int index = 0;
					int next_i = 0;
					for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
					{
						int i2 = *it;
						if(prefix1[i2] > prev_score)
						{
							prev_score = prefix1[i2];
							index = 1;
							next_i = i2;
						}
					}
					int i2 = prev[i];
					if(i2 != -1 && prefix1_L[i2] > prev_score)
					{
						prev_score = prefix1_L[i2];
						index = 2;
						next_i = i2;
					}

					i = next_i;
					if(index == 0)
					{
						break;
					}
					else if (index == 1)
					{
				// t is unchanged
					}
					else
					{
						t = cell_prefix1_L;
					}
				}
				else if(t == cell_suffix1)
				{
					path1.push_back(N-j);
					if (s == 0)
					{
						double next_score = 0.0;
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if(suffix2[N-j2]+ptm_penalty > next_score)
							{
								next_score = suffix2[N-j2]+ptm_penalty;
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = next[N-j];
						if(j2 != -1 && suffix2_R[N-j2]+ptm_penalty > next_score)
						{
							next_score = suffix2_R[N-j2]+ptm_penalty;
							index = 2;
							next_j = N-j2;
						}
						for(std::vector<int>::const_iterator it = right_jumps2[N-j].begin(); it != right_jumps2[N-j].end(); ++it)
						{
							int j2 = *it;
							if(suffix1max[N-j2] > next_score)
							{
								next_score = suffix1max[N-j2];
								index = 3;
								next_j = N-j2;
							}
						}
						j2 = next2[N-j];
						if(j2 != -1 && suffix1_R[N-j2] > next_score)
						{
							next_score = suffix1_R[N-j2];
							index = 4;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							break;
						}
						else if(index == 1)
						{
							t = cell_suffix2;
						}
						else if(index == 2)
						{
							t = cell_suffix2_R;
						}
						else if(index == 3)
						{
							// t is unchanged
							s = (std::max_element(suffix1[j].begin(),suffix1[j].end())-suffix1[j].begin());
						}
						else
						{
							t = cell_suffix1_R;
						}
					}
					else
					{
						// t is unchanged
						j = right_neighbors[N-j][s-1];
						j = N-j;
						s = (std::max_element(suffix1[j].begin(),suffix1[j].end())-suffix1[j].begin());
					}
				}
				else if(t == cell_prefix2)
				{
					path2.insert(path2.begin(), i);
					if(s == 0)
					{
						double prev_score = 0.0;
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(prefix1[i2]+ptm_penalty > prev_score)
							{
								prev_score = prefix1[i2]+ptm_penalty;
								index = 1;
								next_i = i2;
							}
						}
						int i2 = prev[i];
						if(i2 != -1 && prefix1_L[i2]+ptm_penalty > prev_score)
						{
							prev_score = prefix1_L[i2]+ptm_penalty;
							index = 2;
							next_i = i2;
						}
						for(std::vector<int>::const_iterator it = left_jumps2[i].begin(); it != left_jumps2[i].end(); ++it)
						{
							int i2 = *it;
							if(prefix2max[i2] > prev_score)
							{
								prev_score = prefix2max[i2];
								index = 3;
								next_i = i2;
							}
						}
						i2 = prev2[i];
						if(i2 != -1 && prefix2_L[i2] > prev_score)
						{
							prev_score = prefix2_L[i2];
							index = 4;
							next_i = i2;
						}

						i = next_i;
						if(index == 0)
						{
							break;
						}
						else if(index == 1)
						{
							t = cell_prefix1;
						}
						else if(index == 2)
						{
							t = cell_prefix1_L;
						}
						else if(index == 3)
						{
							// t is unchanged
							s = (std::max_element(prefix2[i].begin(),prefix2[i].end())-prefix2[i].begin());
						}
						else
						{
							t = cell_prefix2_L;
						}
					}
					else
					{
						// t is unchanged
						i = left_neighbors[i][s-1];
						s = (std::max_element(prefix2[i].begin(),prefix2[i].end())-prefix2[i].begin());
					}
				}
				else if(t == cell_suffix2)
				{
					path2.push_back(N-j);
					double next_score = 0.0;
					int index = 0;
					int next_j = 0;
					for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
					{
						int j2 = *it;
						if (suffix2[N-j2] > next_score)
						{
							next_score = suffix2[N-j2];
							index = 1;
							next_j = N-j2;
						}
					}
					int j2 = next[N-j];
					if (j2 != -1 && suffix2_R[N-j2] > next_score) {
					next_score = suffix2_R[N-j2];
					index = 2;
					next_j = N-j2;
				}

				j = next_j;
				if(index == 0)
				{
					break;
				}
				else if(index == 1)
				{
					// t is unchanged
				}
				else
				{
					t = cell_suffix2_R;
				}
			}
				else if(t == cell_prefix1_L)
				{
					if(prefix1_L[i] == prefix1[i])
					{
						t = cell_prefix1;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == cell_prefix2_L)
				{
					if(prefix2_L[i] == prefix2max[i])
					{
						t = cell_prefix2;
						s = (std::max_element(prefix2[i].begin(),prefix2[i].end())-prefix2[i].begin());
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == cell_suffix1_R)
				{
					if (suffix1_R[j] == suffix1max[j])
					{
						t = cell_suffix1;
						s = (std::max_element(suffix1[j].begin(),suffix1[j].end())-suffix1[j].begin());
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == cell_suffix2_R)
				{
					if(suffix2_R[j] == suffix2[j])
					{
						t = cell_suffix2;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == cell_D1)
				{
					if(i >= j)
					{
						path1.insert(path1.begin(), i);
						double prev_score = suffix2[j]+ptm_penalty;
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(D1[i2][j] > prev_score)
							{
								prev_score = D1[i2][j];
								index = 1;
								next_i = i2;
							}
						}
						int i2 = prev[i];
						if(i2 != -1 && M1_L[i2][j] > prev_score)
						{
							prev_score = M1_L[i2][j];
							index = 2;
							next_i = i2;
						}

						i = next_i;
						// note that if index=0, then the value of i is irrelevant.
						if(index == 0)
						{
							t = cell_suffix2;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = cell_M1_L;
						}
					}
					else
					{
						path2.push_back(N-j);
						double next_score = prefix1[i]+ptm_penalty;
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if (D1[i][N-j2] > next_score)
							{
								next_score = D1[i][N-j2];
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = next[N-j];
						if(j2 != -1 && M1_R[i][N-j2] > next_score)
						{
							next_score = M1_R[i][N-j2];
							index = 2;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							t = cell_prefix1;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = cell_M1_R;
						}
					}
				}
				else if(t == cell_D2)
				{
					if (i > j)
					{
						path2.insert(path2.begin(), i);
						if(s == 0)
						{
							double prev_score = suffix2[j];
							int index = 0;
							int next_i = 0;
							for (std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
							{
								int i2 = *it;
								if(D1[i2][j] > prev_score)
								{
									prev_score = D1[i2][j];
									index = 1;
									next_i = i2;
								}
							}
							int i2 = prev[i];
							if(i2 != -1 && M1_L[i2][j] > prev_score)
							{
								prev_score = M1_L[i2][j];
								index = 2;
								next_i = i2;
							}
							for(std::vector<int>::const_iterator it = left_jumps2[i].begin(); it != left_jumps2[i].end(); ++it)
							{
								int i2 = *it;
								if(D2max[i2][j] > prev_score)
								{
									prev_score = D2max[i2][j];
									index = 3;
									next_i = i2;
								}
							}
							i2 = prev2[i];
							if(i2 != -1 && M2_L[i2][j] > prev_score)
							{
								prev_score = M2_L[i2][j];
								index = 4;
								next_i = i2;
							}

							i = next_i;
							if(index == 0)
							{
								t = cell_suffix2;
							}
							else if(index == 1)
							{
								t = cell_D1;
							}
							else if(index == 2)
							{
								t = cell_M1_L;
							}
							else if(index == 3)
							{
								// t is unchanged
								s = (std::max_element(D2[i][j].begin(),D2[i][j].end())-D2[i][j].begin());
							}
							else
							{
								t = cell_M2_L;
							}
						}
						else
						{
							// t is unchanged
							i = left_neighbors[i][s-1];
							s = (std::max_element(D2[i][j].begin(),D2[i][j].end())-D2[i][j].begin());
						}
					}
					else
					{ // i <= j
						path2.push_back(N-j);
						double next_score = prefix2[i][s];
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if(D2[i][N-j2][s] > next_score)
							{
								next_score = D2[i][N-j2][s];
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = next[N-j];
						if(j2 != -1 && M2_R[i][N-j2][s] > next_score)
						{
							index = 2;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							t = cell_prefix2;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = cell_M2_R;
						}
					}
				}
				else if(t == cell_D3)
				{
					if(i >= j)
					{
						path1.insert(path1.begin(), i);
						double prev_score = suffix1[j][s];
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = left_jumps[i].begin(); it != left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(D3[i2][j][s] > prev_score)
							{
								prev_score = D3[i2][j][s];
								index = 1;
								next_i = i2;
							}
						}
						int i2 = prev[i];
						if(i2 != -1 && M3_L[i2][j][s] > prev_score)
						{
							index = 2;
							next_i = i2;
						}

						i = next_i;
						if(index == 0)
						{
							t = cell_suffix1;
						}
						else if (index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = cell_M3_L;
						}
					}
					else
					{ // i < j
						path1.push_back(N-j);
						if(s == 0)
						{
							double next_score = prefix1[i];
							int index = 0;
							int next_j = 0;
							for(std::vector<int>::const_iterator it = right_jumps[N-j].begin(); it != right_jumps[N-j].end(); ++it)
							{
								int j2 = *it;
								if(D1[i][N-j2] > next_score)
								{
									next_score = D1[i][N-j2];
									index = 1;
									next_j = N-j2;
								}
							}
							int j2 = next[N-j];
							if(j2 != -1 && M1_R[i][N-j2] > next_score)
							{
								next_score = M1_R[i][N-j2];
								index = 2;
								next_j = N-j2;
							}
							for(std::vector<int>::const_iterator it = right_jumps2[N-j].begin(); it != right_jumps2[N-j].end(); ++it)
							{
								int j2 = *it;
								if(D3max[i][N-j2] > next_score)
								{
									next_score = D3max[i][N-j2];
									index = 3;
									next_j = N-j2;
								}
							}
							j2 = next2[N-j];
							if(j2 != -1 && M3_R[i][N-j2] > next_score)
							{
								next_score = M3_R[i][N-j2];
								index = 4;
								next_j = N-j2;
							}

							j = next_j;
							if(index == 0)
							{
								t = cell_prefix1;
							}
							else if(index == 1)
							{
								t = cell_D1;
							}
							else if(index == 2)
							{
								t = cell_M1_R;
							}
							else if( index == 3)
							{
								// t is unchanged
								s = (std::max_element(D3[i][j].begin(),D3[i][j].end())-D3[i][j].begin());
							}
							else
							{
								t = cell_M3_R;
							}
						}
						else
						{
							// t is unchanged
							j = right_neighbors[N-j][s-1];
							j = N-j;
							s = (std::max_element(D3[i][j].begin(),D3[i][j].end())-D3[i][j].begin());
						}
					}
				}
				else if(t == cell_M1_L)
				{
					if(M1_L[i][j] == D1[i][j])
					{
						t = cell_D1;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == cell_M2_L)
				{
					if (M2_L[i][j] == D2max[i][j])
					{
						t = cell_D2;
						s = (std::max_element(D2[i][j].begin(),D2[i][j].end())-D2[i][j].begin());
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == cell_M3_L)
				{
					if (M3_L[i][j][s] == D3[i][j][s])
					{
						t = cell_D3;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == cell_M1_R)
				{
					if(M1_R[i][j] == D1[i][j])
					{
						t = cell_D1;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == cell_M2_R)
				{
					if(M2_R[i][j][s] == D2[i][j][s])
					{
						t = cell_D2;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == cell_M3_R)
				{
					if(M3_R[i][j] == D3max[i][j])
					{
						t = cell_D3;
						s = (std::max_element(D3[i][j].begin(),D3[i][j].end())-D3[i][j].begin());
					}
					else
					{
						j -= 1;
					}
				}
			}
			return std::pair<double, std::pair<std::vector<int>,std::vector<int> > >(best_score/parent_mass1, std::pair<std::vector<int>,std::vector<int> >(path1, path2));
		}
		///

	};
	///

	/**
		@brief Aligns the peaks of two spectra

		@htmlinclude OpenMS_AntisymetricAlignment.parameters

		@ingroup SpectraComparison

	/// @todo debug output in macro
	*/
	template <typename PeakT = Peak1D>
	class AntisymetricAlignment
		: public DefaultParamHandler
	{
		typedef MSSpectrum<PeakT> SpectrumType;

	public:

		// @name Constructors and Destructors
		// @{
		/// default constructor
		AntisymetricAlignment()
		: DefaultParamHandler("AntisymetricAlignment")
		{
				defaults_.setValue("peak_tolerance", 0.3, "Defines the absolut (in Da) peak tolerance");
				defaults_.setValue("min_dist", 57.0214637230 , "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
				/// @improvement DP penalty scores here
				defaultsToParam_();
		}

		/// copy constructor
		AntisymetricAlignment(const AntisymetricAlignment& source)
			: DefaultParamHandler(source)
		{
		}

		/// destructor
		virtual ~AntisymetricAlignment()
		{
		}

		/// assignment operator
		AntisymetricAlignment& operator = (const AntisymetricAlignment& source)
		{
			if (this != &source)
			{
				DefaultParamHandler::operator = (source);
			}
			return *this;
		}

		// @}

		/**
			@brief Will calculate valid match jumps left and right to a given peak in a given spectrum

			@param index to given peak
			@param s1sym the given spectrum
			@param jump_masses the valid mass jumps for the matches
			@param left_matches will contain the valid match jumps to the left of index
			@param right_matches will contain the valid match jumps to the right of index
		*/
		void findMatchingJumps(Size& index, SpectrumType& s1sym, std::vector<Real>& jump_masses,
						std::vector<int>& leftMatches, std::vector<int>& rightMatches) const
		{
			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");
			Real pt_plus = peak_tolerance +.00001f;
			Size match_counter = 0;
			DoubleReal current_diff;
			DoubleReal max_jump = jump_masses[jump_masses.size()-1]+peak_tolerance;

			//~ if we start leftmost this can be skipped as there is nothing to jump left to
			if(index>0)
			{
				//~ max |index| left jumps possible
				leftMatches.resize(index);
				//~ but so far none
				match_counter=0;

				//~ typename SpectrumType::Iterator it_lo1 = s1sym.begin()+index-1;
				//~ typename SpectrumType::Iterator it_lo2 = s1sym.MZBegin(s1sym[index].getMZ()-max_jump-.00001);

				std::reverse_iterator< typename SpectrumType::Iterator > it_lo1(s1sym.begin()+index); //points to before 'index'
				std::reverse_iterator< typename SpectrumType::Iterator > it_lo2(s1sym.MZBegin(s1sym[index].getMZ()-max_jump-.00001)); //points to before MZBegin... maybe before s1sym.begin() -std::min(s1sym[0].getMZ(), MZBegin...)?
				std::vector<Real>::iterator jm_left_pivot = jump_masses.begin();

				/// @improvement posssible speedup or slowdown by comp if it_lo2.getMZ() <= it_lo1.getMZ() and only go on then ...
				/*debug std::cout << " it_lo1: " << it_lo1->getMZ() << " reverse until incl. :  " << (it_lo2-1)->getMZ() << std::endl ; */

				//~ iterate from lo1 to lo2(not included!) and find if mz diff between lo1(iterated) and index concurs with a entry in jump_masses +-peak_tolerance
				for(;it_lo1 != it_lo2; ++it_lo1)
				{
					current_diff = s1sym[index].getMZ()-it_lo1->getMZ();
					/*debug  std::cout << " current diff: " << current_diff; */
					//~ difference will only grow so we can advance the pivot
					jm_left_pivot = std::lower_bound(jm_left_pivot, jump_masses.end(), (Real)current_diff-pt_plus);

					if(jm_left_pivot!=jump_masses.end() and EqualInTolerance<Real>(pt_plus)((Real)current_diff,*jm_left_pivot))
					{
						leftMatches[match_counter++]=(s1sym.size()-1)-(it_lo1-s1sym.rbegin());
						/*debug  std::cout << " -left jo " ;*/
					}
				}
				/*debug  std::cout << '\n';*/
				leftMatches.resize(match_counter);
			}

			if(index<(s1sym.size()-1))
			{
				//~ max |s1sym|-|index| right jumps possible
				rightMatches.resize(s1sym.size()-index);
				//~ but so far none
				match_counter=0;

				typename SpectrumType::Iterator it_hi1 = s1sym.begin()+index+1;
				typename SpectrumType::Iterator it_hi2 = s1sym.MZEnd(s1sym[index].getMZ()+max_jump+.00001);
				std::vector<Real>::iterator jm_right_pivot = jump_masses.begin();
				/*debug std::cout << " it_hi1: " << it_hi1->getMZ() << " until incl. :  " << (it_hi2-1)->getMZ() << std::endl ; */

				//~ iterate from hi1 to <hi2 and find if mz diff between it_hi1(iterated) and index concurs with a entry in jump_masses +-peak_tolerance
				for(;it_hi1 != it_hi2; ++it_hi1)
				{
					current_diff = it_hi1->getMZ()-s1sym[index].getMZ();
					/*debug  std::cout << " current diff: " << current_diff;*/
					//~ difference will only grow so we can advance the pivot
					jm_right_pivot = std::lower_bound(jm_right_pivot, jump_masses.end(), (Real)current_diff-pt_plus);

					if(jm_right_pivot!=jump_masses.end() and EqualInTolerance<Real>(pt_plus)((Real)current_diff,*jm_right_pivot))
					{
						rightMatches[match_counter++]=it_hi1-s1sym.begin();
						/*debug  std::cout << " -right jo " ;*/
					}
				}
				/*debug  std::cout << '\n';*/
				rightMatches.resize(match_counter);
			}
		}
		///

		/**
			@brief Will yield a symetric spectrum

			@param spectrum to be made symetric

			@return a symetric version of spectrum with a integer MetaDataArray, where x[i]>0 indicates a synthetic peak

			result is spectrum with property sym[i].getMZ + sym[n-i-1].getMZ = parentMass-Constants::PROTON_MASS_U, adding entries setIntensity zero if necessary
		@important previous sortByPosition() of spectrum
		@important is subtraction of a H enough?? how to decide between prm and exp!?
		@important is the 0 and pm peak addition neccessary!?
		*/
		MSSpectrum<PeakT> getSymetricSpectrum(const MSSpectrum<PeakT>& spectrum) const
		{
			SpectrumType sym(spectrum);
			if(sym.empty())
			{
				return sym;
			}

			typename SpectrumType::IntegerDataArray ida_symetric; ida_symetric.resize(sym.size(),0);
			sym.getIntegerDataArrays().push_back(ida_symetric);
			sym.getIntegerDataArrays().push_back(ida_symetric);

			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");

			DoubleReal pm = sym.getPrecursors().front().getMZ();
			int pm_charge = sym.getPrecursors().front().getCharge();
			/// @important singly charged mass difference!
			pm = pm * pm_charge - (pm_charge-1)*Constants::PROTON_MASS_U;
			/* debug std::cout << "Sym pm "<< pm << std::endl; */


			Size sym_pasttheend_index = sym.size();
			for(Size sym_index = 0; sym_index < sym_pasttheend_index; ++sym_index)
			{
				if(sym.getIntegerDataArrays().front()[sym_index]!=0)
				{
					continue;
				}
				DoubleReal target = pm-sym[sym_index].getMZ()+Constants::PROTON_MASS_U;

				typename SpectrumType::Iterator lo = sym.MZBegin(target-peak_tolerance);
				typename SpectrumType::Iterator hi = sym.MZEnd(lo,target+peak_tolerance,sym.begin()+sym_pasttheend_index);
				/*debug std::cout << target-peak_tolerance << " - " << target+peak_tolerance << std::endl; */

				//~ pairs get metadataarray entries -n/n where n is number of the pair (count from 1) and the lower mz peak has the negative pair number
				int sign_it_mz = -1;
				//~ int sign_target = sign_it_mz * -1
				((target-sym[sym_index].getMZ())<0)?sign_it_mz*=(-1):sign_it_mz=sign_it_mz;

				//~ empty range - add a peak
				if(lo==hi)
				{
					/// @improvement catch target-is-negative case
					PeakT tmp;
					tmp.setPosition(target);
					sym.push_back(tmp);
					sym.getIntegerDataArrays().front().push_back( sign_it_mz*(-1)* sym_index+1 );
					sym.getIntegerDataArrays().front()[sym_index]= sign_it_mz * sym_index+1;
					sym.getIntegerDataArrays()[1].push_back(1);
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::Iterator it_mirror_mz = lo; it_mirror_mz != hi; ++it_mirror_mz)
					{
						if( fabs(it_mirror_mz->getMZ()-sym[sym_index].getMZ()) < fabs(lo->getMZ()-sym[sym_index].getMZ()) )
						{
							lo=it_mirror_mz;
						}
					}
					sym.getIntegerDataArrays().front()[sym_index]= sign_it_mz*(sym_index+1);
					sym.getIntegerDataArrays().front()[(lo-sym.begin())]= sign_it_mz*(-1)*(sym_index+1);
				}
			}

			std::swap(sym.getIntegerDataArrays()[0],sym.getIntegerDataArrays()[1]);
			sym.getIntegerDataArrays()[0].setName("synthetic peaks");
			sym.getIntegerDataArrays().resize(1);
			sym.sortByPosition();

			return sym;
		}
		///

		/**
			@brief Method to calculate a antisymetric alignment

			@param res_1 Spectrum containing the aligned peaks from s1
			@param res_2 Spectrum containing the aligned peaks from s2
			@param score containing the score (best_score/parent_mass s1)
			@param mod_pos the position in res_1 that has to be 'modified' (even if it is the aligned counterpart in res2 that had a subtraction)
			@param s2 the first input spectrum
			@param s2 the second input spectrum

			This method must have the first spectrum have a lesser or equal parent mass than the second spectrum, else a error is thrown. Also the first spectrum must be symetric, so synthetic peaks with zero intensity might be inserted. Input does not. Also the resulting spectra are intended for consensus making which takes a list of spectra.
		*/
		void getAntisymetricAlignment(MSSpectrum<PeakT>& res_1, MSSpectrum<PeakT>& res_2, DoubleReal& score, DoubleReal& mod_pos, const MSSpectrum<PeakT>& s1, const MSSpectrum<PeakT>& s2) const
		{
			//~ float sameVertexPenalty=-5, float ptmPenalty=-5, bool forceSymmetry=true, bool addZPMmatches=false

			if (s1.getPrecursors().front().getMZ() > s2.getPrecursors().front().getMZ())
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "prerequisite is the s1-precursor mz is not greater than s2-precursor mz");
			}

			//~ res_1.clear(true);
			//~ res_2.clear(true);
			res_1.clear(); res_1.getFloatDataArrays().clear();res_1.getIntegerDataArrays().clear();res_1.getStringDataArrays().clear();
			res_2.clear(); res_2.getFloatDataArrays().clear();res_2.getIntegerDataArrays().clear();res_2.getStringDataArrays().clear();


			if(s1.empty() or s2.empty())
			{
				return;
			}

			//~ attention: max. number of AA jumped is always 1

			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");
			//~ Real parentmass_tolerance = (Real)param_.getValue("parentmass_tolerance");
			//~ Real max_shift = (Real)param_.getValue("max_shift");
			Real min_dist = (Real)param_.getValue("min_dist");
			DoubleReal pm_s1 = s1.getPrecursors().front().getMZ();
			int c_s1 = s1.getPrecursors().front().getCharge();
			DoubleReal pm_s2 = s2.getPrecursors().front().getMZ();
			int c_s2 = s2.getPrecursors().front().getCharge();
			/// @important singly charged mass difference!
			Real pm_diff = (pm_s2*c_s2 + (c_s2-1)*Constants::PROTON_MASS_U)-(pm_s1*c_s1 + (c_s1-1)*Constants::PROTON_MASS_U);

			//jump_masses[i] holds internal mass of aa i
			//sort jump_masses
			//make unique in resolution tolerance (e.g. cause of Q & K)
			const ResidueDB* res_db = ResidueDB::getInstance();
			std::vector<Real> jump_masses;
			for(ResidueDB::ResidueConstIterator it = res_db->beginResidue(); it != res_db->endResidue(); ++it)
			{
				DoubleReal w = (*it)->getMonoWeight() - (*it)->getInternalToFullMonoWeight();
				/// @important: this makes only singly charged ion peaks to potential jump matches and therefor potential alignment-partners (so a test is only reasonable with singly charged ion peaks in da spec)
				if(w>0)
				{
					jump_masses.push_back(w);
				}
			}

			std::sort(jump_masses.begin(),jump_masses.end());
			std::vector<Real>::iterator end_new = std::unique(jump_masses.begin(),jump_masses.end(), EqualInTolerance<Real>(peak_tolerance));
			jump_masses.resize( end_new - jump_masses.begin() );
			DoubleReal jumps_supremum=0;
			jumps_supremum = jump_masses.back()+2*peak_tolerance+0.00001;
			for(Size w = 0; w < jump_masses.size(); ++w)
			{
				/* debug  std::cout << jump_masses[w] << ", ";*/
			}
			/* debug  std::cout << '\n';*/

			//~ change s1 to align antisymetrical into res_1
			res_1 = getSymetricSpectrum(s1);
			//~ indices in s2 for normal matches
			/* debug  std::cout << "symetric size is: " << res_1.size() << std::endl;*/

			typename SpectrumType::IntegerDataArray ida_matches; ida_matches.resize(res_1.size(),-1);
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin()+1, ida_matches);
			res_1.getIntegerDataArrays()[1].setName("normal match indices");
			//~ indices in s2 for pm_diff shift matches
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin()+2, ida_matches);
			res_1.getIntegerDataArrays()[2].setName("shift match indices");

			//~ calculate common/common2 (StringDataArrays) and common_scores/common2_scores (intensities)
			for(typename SpectrumType::Iterator it_res = res_1.begin(); it_res != res_1.end(); ++it_res)
			{
				bool no_match_normal = false;
				typename SpectrumType::ConstIterator lo = s2.MZBegin(it_res->getMZ()-peak_tolerance);
				typename SpectrumType::ConstIterator hi = s2.MZEnd(it_res->getMZ()+peak_tolerance);
				if(lo==hi)
				{
					no_match_normal = true;
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::ConstIterator it_match_mz = lo; it_match_mz != hi; ++it_match_mz)
					{
						if( fabs(it_match_mz->getMZ()-it_res->getMZ()) < fabs(lo->getMZ()-it_res->getMZ()) )
						{
							lo=it_match_mz;
						}
					}
					res_1.getIntegerDataArrays()[1][(it_res-res_1.begin())] = lo-s2.begin();
				}

				lo = s2.MZBegin(it_res->getMZ()+pm_diff-peak_tolerance);
				hi = s2.MZEnd(it_res->getMZ()+pm_diff+peak_tolerance);
				if(lo==hi)
				{
					if(no_match_normal)
					{
						it_res->setIntensity(0);
					}
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::ConstIterator it_match_mz = lo; it_match_mz != hi; ++it_match_mz)
					{
						if( fabs(it_match_mz->getMZ()-it_res->getMZ()) < fabs(lo->getMZ()-it_res->getMZ()) )
						{
							lo=it_match_mz;
						}
					}
					res_1.getIntegerDataArrays()[2][(it_res-res_1.begin())] = lo-s2.begin();
				}
			}

			/* debug
			std::cout << "§ " ;
			for(Size i = 0; i < res_1.getIntegerDataArrays()[1].size(); ++i)
			{
				std::cout << res_1.getIntegerDataArrays()[1][i];
			}
			std::cout << "$ ";
			for(Size i = 0; i < res_1.getIntegerDataArrays()[2].size(); ++i)
			{
				std::cout << res_1.getIntegerDataArrays()[2][i];
			}
			std::cout << std::endl;
			*/

			//~ number of complementing peaks in res_1 that both do not match the right peaks in s2
			Size removed_pairs = 0;
			//~ remove the peak pairs mentioned above
			Size res_1_index = 0;
			while((Size)res_1_index < (Size)((res_1.size()/2) - removed_pairs) )
			{
				if(res_1[res_1_index].getIntensity() == 0)
				{
					Size res_mirror_index = res_1.size() - 1 - res_1_index;
					if(res_1[res_mirror_index].getIntensity() == 0 )
					{
						for(Size i = 0; i < 3; ++i)
						{
							typename SpectrumType::IntegerDataArray::iterator it_mirror_ida = res_1.getIntegerDataArrays()[i].begin()+res_mirror_index;
							res_1.getIntegerDataArrays()[i].erase(it_mirror_ida);
						}
						for(Size i = 0; i < 3; ++i)
						{
							typename SpectrumType::IntegerDataArray::iterator it_ida = res_1.getIntegerDataArrays().front().begin()+res_1_index;
							res_1.getIntegerDataArrays()[i].erase(it_ida);
						}
						typename SpectrumType::Iterator res_mirror_it = res_1.begin() + res_mirror_index;
						res_1.erase(res_mirror_it);
						typename SpectrumType::Iterator res_1_it = res_1.begin() + res_1_index;
						res_1.erase(res_1_it);
						++removed_pairs;
					}
					else
					{
						++res_1_index;
					}
				}
				else
				{
					++res_1_index;
				}
			}

			std::vector<DoubleReal> peaks(res_1.size());  for(Size i=0;i<res_1.size();++i) peaks[i]=res_1[i].getMZ();
			std::vector<DoubleReal> peaks2(s2.size());   for(Size i=0;i<s2.size();++i) peaks2[i]=s2[i].getMZ();
			std::vector<int> common(res_1.size(),-1);
			std::vector<int> common2(res_1.size(),-1);
			std::vector<double> common_scores(res_1.size(),0);
			std::vector<double> common2_scores(res_1.size(),0);

			//~ fill the above from res
			for(Size i=0;i<res_1.size();++i)
			{
				if(res_1.getIntegerDataArrays()[1][i]>-1)
				{
					common[i] = res_1.getIntegerDataArrays()[1][i];
					/* debug  std::cout << "common "<< i << " = " << common[i] << std::endl;*/
					common_scores[i] = res_1[i].getIntensity()+s2[common[i]].getIntensity();
					/* debug  std::cout << "common_scores "<< i << " = " << common_scores[i] << std::endl;*/
				}
				else if(res_1.getIntegerDataArrays()[2][i]>-1)
				{
					common2[i] = res_1.getIntegerDataArrays()[2][i];
					/* debug  std::cout << "common2 "<< i << " = " << common2[i] << std::endl;*/
					common2_scores[i] = res_1[i].getIntensity()+s2[common2[i]].getIntensity();
					/* debug  std::cout << "common2_scores "<< i << " = " << common2_scores[i] << std::endl;*/
				}
			}

			std::vector<int> prev;      prev.resize(res_1.size(), -1);
			std::vector<int> next;      next.resize(res_1.size(),-1);
			std::vector<int> prev2;     prev2.resize(res_1.size(),-1);
			std::vector<int> next2;     next2.resize(res_1.size(),-1);
			std::vector<std::vector<int> > left_jumps(res_1.size());
			std::vector<std::vector<int> > left_jumps2(res_1.size());
			std::vector<std::vector<int> > right_jumps(res_1.size());
			std::vector<std::vector<int> > right_jumps2(res_1.size());
			std::vector<std::vector<int> > left_neighbors(res_1.size());
			std::vector<std::vector<int> > right_neighbors(res_1.size());

			float delta = std::min(std::max(min_dist-peak_tolerance,pm_diff+peak_tolerance),min_dist-peak_tolerance*2.0f);
			float delta2 = jumps_supremum-peak_tolerance;

			/* debug  std::cout << "delta "<< delta << " delta2 " << delta2 << std::endl;*/

			//~ fill the above from res_1/peaks
			/* debug  std::cout << " * "<< res_1.size() << std::endl;*/
			for(Size i = 0; i < res_1.size(); ++i)
			{
				/* debug  std::cout << " i: "<< i << std::endl;*/
				DoubleReal lo1 = res_1[i].getMZ()-delta2;
				if(lo1<=0)
				{
					prev[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_lo1 = res_1.MZEnd(res_1.begin(),lo1,res_1.begin()+i);
					it_lo1!=res_1.begin()?prev[i]=((it_lo1-res_1.begin())-1):prev[i]=-1;
				}
				DoubleReal lo2 = res_1[i].getMZ()-std::max(delta,delta2);
				if(lo2<=0)
				{
					prev2[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_lo2 = res_1.MZEnd(lo2);
					it_lo2!=res_1.begin()?prev2[i]=((it_lo2-res_1.begin())-1):prev2[i]=-1;
				}

				DoubleReal hi1 = res_1[i].getMZ()+ delta2;
				if(hi1>=res_1.rbegin()->getMZ())
				{
					next[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_hi1 = res_1.MZBegin(res_1.begin()+i,hi1,res_1.end());
					it_hi1!=res_1.end()?next[i]=(it_hi1-res_1.begin()):next[i]=-1;
				}
				DoubleReal hi2 = res_1[i].getMZ()+std::max(delta,delta2);
				if(hi2>=res_1.rbegin()->getMZ())
				{
					next2[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_hi2 = res_1.MZBegin(res_1.begin()+i,hi2,res_1.end());
					it_hi2!=res_1.end()?next2[i]=(it_hi2-res_1.begin()):next2[i]=-1;
				}

				/* debug  std::cout << " prev: "<< prev[i] << " prev2: " << prev2[i] << std::endl;*/
				/* debug  std::cout << " next: "<< next[i] << " next2: " << next2[i] << std::endl;*/

				findMatchingJumps(i, res_1, jump_masses, left_jumps[i], right_jumps[i]);
				Size neighbor_count=0, jumps_count=0;
				left_neighbors[i].resize(left_jumps[i].size());
				left_jumps2[i].resize(left_jumps[i].size());
				/* debug  std::cout << "left jumps size "<< left_jumps[i].size() << " ";*/
				for(Size j=0; j<left_jumps[i].size(); ++j)
				{
					if(peaks[i]-peaks[left_jumps[i][j]]<=delta)
					{
						left_neighbors[i][neighbor_count++]=left_jumps[i][j];
						/* debug  std::cout << " ln] "<< left_jumps[i][j] << " ";*/
					}
					if(peaks[i]-peaks[left_jumps[i][j]]>=delta)
					{
						left_jumps2[i][jumps_count++]=left_jumps[i][j];
						/* debug  std::cout << " lj] "<< left_jumps[i][j] << " ";*/
					}
				}
				left_neighbors[i].resize(neighbor_count);
				left_jumps2[i].resize(jumps_count);
				/* debug */ std::cout << std::endl;

				neighbor_count=0, jumps_count=0;
				right_neighbors[i].resize(right_jumps[i].size());
				right_jumps2[i].resize(right_jumps[i].size());
				/* debug  std::cout << "right jumps size "<< right_jumps[i].size() << " ";*/
				for(Size j=0; j<right_jumps[i].size(); ++j)
				{
					if(peaks[right_jumps[i][j]]-peaks[i]<=delta)
					{
						right_neighbors[i][neighbor_count++]=right_jumps[i][j];
						/* debug  std::cout << " rn] "<< right_jumps[i][j] << " ";*/
					}
					if(peaks[right_jumps[i][j]]-peaks[i]>=delta)
					{
						right_jumps2[i][jumps_count++]=right_jumps[i][j];
						/* debug  std::cout << " rj] "<< right_jumps[i][j] << " ";*/
					}
				}
				right_neighbors[i].resize(neighbor_count);
				right_jumps2[i].resize(jumps_count);
				/* debug  std::cout << std::endl;*/
			}

			res_1.clear();
			res_2 = s2;
			res_2.clear();
			//~ implicitly copied (and not erased):
			//~ Precursors,RT,MSLevel, MetaDataArrays

			/* debug  std::cout << " aligning" << std::endl;*/

			/// @improvement btw is this right with massMH = H2O + Constants::PROTON_MASS_U?!?!?!
			//~ static const double massMH = EmpiricalFormula("H2O").getMonoWeight() + Constants::PROTON_MASS_U;
			/// @improvement find a good penalty adjustment
			DoubleReal sameVertexPenalty = -5, ptmPenalty = -5;
			//~ DoubleReal sameVertexPenalty = -1000000, ptmPenalty = -200;
			//~ DoubleReal sameVertexPenalty = 0, ptmPenalty = 0;

			std::pair<double, std::pair< std::vector<int>,std::vector<int> > > asym_align;

			AntisymetricDP dp;
			asym_align = dp.calculateAlignementPath(/* pm_s1 - massMH */ (pm_s1*c_s1 + (c_s1-1)*Constants::PROTON_MASS_U), /* pm_s2 - massMH */(pm_s2*c_s2 + (c_s2-1)*Constants::PROTON_MASS_U), peaks, peaks2, common, common2, common_scores, common2_scores, prev, next, prev2, next2, left_jumps, right_jumps, left_jumps2, right_jumps2, left_neighbors, right_neighbors, peak_tolerance, sameVertexPenalty, ptmPenalty);

			score = asym_align.first;

			Size num_aligned_peaks = asym_align.second.first.size()+asym_align.second.second.size();

			/* debug  std::cout << " aligned: "<< asym_align.second.first.size() << " + " << asym_align.second.second.size() << std::endl;*/

			typename SpectrumType::IntegerDataArray ida_synthetic; ida_synthetic.resize(num_aligned_peaks,0);
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin(),ida_synthetic);
			(res_1.getIntegerDataArrays().begin())->setName("synthetic peaks");
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin()+1,ida_synthetic);
			(res_1.getIntegerDataArrays().begin()+1)->setName("modification position");
			res_2.getIntegerDataArrays().insert(res_2.getIntegerDataArrays().begin(), ida_synthetic);
			(res_2.getIntegerDataArrays().begin())->setName("synthetic peaks");
			res_2.getIntegerDataArrays().insert(res_2.getIntegerDataArrays().begin()+1, ida_synthetic);
			(res_2.getIntegerDataArrays().begin()+1)->setName("modification position");

			/* debug  std::cout << " determine mod_pos: ";*/
			if(asym_align.second.first.size()==0)
			{
				mod_pos=0;  // No peaks in idx1 => mod at the start
				(res_1.getIntegerDataArrays().begin()+1)->at(0)=-1;
				(res_2.getIntegerDataArrays().begin()+1)->at(0)=-1;
			}
			else
			{
				if(asym_align.second.second.size()==0)
				{
					mod_pos = pm_s1;  // No peaks in idx2 => mod at the end
					(res_1.getIntegerDataArrays().begin()+1)->at(num_aligned_peaks-1)=-1;
					(res_2.getIntegerDataArrays().begin()+1)->at(num_aligned_peaks-1)=-1;
				}
				else
				{
					mod_pos=peaks[asym_align.second.second[0]];  // Otherwise mod was placed at the first mass of res
					(res_1.getIntegerDataArrays().begin()+1)->at(asym_align.second.first.size())=1;
					(res_2.getIntegerDataArrays().begin()+1)->at(asym_align.second.first.size())=1;
				}
			}
			/* debug  std::cout << mod_pos << std::endl;*/

			/* debug  std::cout << " creating spectra: ";*/
			typename SpectrumType::IntegerDataArray& ida_ref = *(res_1.getIntegerDataArrays().begin()+2);
			for(Size i=0; i<asym_align.second.first.size(); ++i)
			{
				PeakT tmp_1, tmp_2;
				tmp_1.setMZ( peaks[asym_align.second.first[i]] );
				if(common[asym_align.second.first[i]]>=0)
				{
					tmp_2.setMZ(peaks2[common[asym_align.second.first[i]]]);
					tmp_2.setIntensity(s2[common[asym_align.second.first[i]]].getIntensity());
					res_2.getIntegerDataArrays().front().push_back(0);
				}
				else
				{
					tmp_2.setMZ(tmp_1.getMZ());
					//~ intensity stays default constructed 0
					res_2.getIntegerDataArrays().front().push_back(1);
				}
				tmp_1.setIntensity(common_scores[asym_align.second.first[i]] - tmp_2.getIntensity());
				(ida_ref[asym_align.second.first[i]]>0)?res_2.getIntegerDataArrays().front()[i]=1:res_2.getIntegerDataArrays().front()[i]=0;
				res_1.push_back(tmp_1);
				res_2.push_back(tmp_2);
			}

			for(Size i=0; i<asym_align.second.second.size(); ++i)
			{
				PeakT tmp_1, tmp_2;
				tmp_1.setMZ( peaks[asym_align.second.second[i]] );
				if(common2[asym_align.second.second[i]]>=0)
				{
					tmp_2.setMZ( peaks2[common2[asym_align.second.second[i]]] );
					tmp_2.setIntensity( s2[common2[asym_align.second.second[i]]].getIntensity() );
					res_2.getIntegerDataArrays().front().push_back(0);
				}
				else
				{
					tmp_2.setMZ(tmp_1.getMZ() + pm_diff);
					//~ intensity stays default constructed 0
					res_2.getIntegerDataArrays().front().push_back(1);
				}
				tmp_1.setIntensity( common2_scores[asym_align.second.second[i]] - tmp_2.getIntensity());
				(ida_ref[asym_align.second.second[i]]>0)?res_2.getIntegerDataArrays().front()[asym_align.second.first.size()+i]=1:res_2.getIntegerDataArrays().front()[asym_align.second.first.size()+i]=0;
				res_1.push_back(tmp_1);
				res_2.push_back(tmp_2);
			}
		}
		///
	};


}
#endif //OPENMS_COMPARISON_SPECTRA_ANTISYMETRICALIGNMENT_H