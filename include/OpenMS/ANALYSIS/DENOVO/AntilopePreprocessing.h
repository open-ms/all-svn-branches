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
// $Maintainer: Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef IDSETUP_H_
#define IDSETUP_H_


#include <vector>
#include <map>
#include <fstream>
#include <cmath>

//OpenMS includes
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>



namespace OpenMS
{

class IdSetup{

	public:

      //TODO  remove interpretation from this class and put it into the scoring function base class to be developed
      struct interpretation
      {
          Int peak_nr;
          Int type_id;
          interpretation(Int a=0,Int b=0):peak_nr(a), type_id(b) {};
      };


    	///@name typedefs
    	//@{      
    	typedef std::vector<UInt>UIntVec;
    	typedef std::vector<DoubleReal>DRealVec;
    	typedef std::vector<String>StringVec;
      typedef std::vector<AASequence>AASeqVec;
    	typedef std::vector<bool>BoolVec;



    	//maps for each integer mass all AA (AA-combinations) matching that mass
    	typedef std::map<UInt,std::vector<String> > AAmap;
    	typedef std::map<UInt,std::vector<String> > StringVecMap;
      typedef std::map<UInt,std::vector<AASequence> > AASeqVecMap;
    	typedef std::map<UInt, String> StringMap;

    	//@}

    	///Constructor
      IdSetup(DoubleReal delta_in, DoubleReal in_precision, bool join_isobarics, bool join_Q_K, UInt max_span ):delta_(delta_in),
                precision_(in_precision),
                join_isobarics_(join_isobarics),
                join_Q_K_(join_Q_K),
                max_span_(max_span)
    	{            
    	};

    	/// default Constructor
    	IdSetup():delta_(0.5),
                precision_(0.1),
                join_isobarics_(true),
                join_Q_K_(true),
                max_span_(3)
    	{       
    	};

      /**method to create array A which is used for generating the edge set --> see chen paper pp.392,
       * the Array is written to a file since it is not problem specific, hence it is not necessary to
       * recreate it everytime  */
      int create_vector_A(BoolVec&, UIntVec&, AASeqVecMap&);

      //void test_io_functions(); // @TODO ist nur zum testen, fliegt raus!

      ///reads the Array A or A_len (type UInt and Bool) from binary files
      template<typename T>
      static std::streampos readVectorFromBinFile(String file, std::streampos pos, std::vector<T> &arr, Int max_elem=-1);

      ///writes the Array A or A_len (type UInt and Bool) to binary files
      template<typename T>
      static std::streampos writeVectorToBinFile(String file, std::streampos pos, const std::vector<T> &arr, Int max_elem=-1);

      static void filterSpectrumInspect(PeakSpectrum &spec)
      {
        Param filter_params;
        filter_params.setValue("windowsize", 56.0);
        filter_params.setValue("peakcount", 3);
        WindowMower win_mower;
        win_mower.setParameters(filter_params);
        win_mower.filterPeakSpectrum(spec);
      };

      ///getters Setter
      DoubleReal getDelta()
      {
        return delta_;
      };

      void setDelta(DoubleReal delta)
      {
        delta_=delta;
      };

      DoubleReal getPrecision()
      {
        return precision_;
      };

      void setPrecision(DoubleReal prec)
      {
        precision_ = prec;
      };

      void setSpan(UInt span)
      {
        max_span_=span;
      }

      UInt getSpan()
      {
        return max_span_;
      }

      bool getIsobaricsFlag()
      {
        return join_isobarics_;
      };

      void setIsobaricsFlag(bool iso)
      {
        join_isobarics_=iso;
      };

      bool getQ_K_Flag()
      {
        return join_Q_K_;
      };

      void setQ_K_Flag(bool q_k_flag)
      {
        join_Q_K_=q_k_flag;
      };

      enum ion_type
      {
        FIRST_TYPE,
        BIon = FIRST_TYPE,
        BIon2,
        BIon_h2o,
        BIon_nh3,
        BIon_h2o_h2o,
        BIon_nh3_h2o,
        YIon,
        YIon2,
        YIon_h2o,
        YIon_nh3,
        YIon_h2o_h2o,
        YIon_nh3_h2o,
        AIon,
        AIon_nh3,
        AIon_h2o,
        CIon,
        ZIon,
        LAST_TYPE = ZIon,
        Invalid
      };

      const static DoubleReal proton_m;

      static const DoubleReal h2o_mass;

      const static DoubleReal nh3_mass;
      const static DoubleReal nh3_h2o_mass;

      const static DoubleReal BIon_offset;
      const static DoubleReal BIon_h2o_offset;
      const static DoubleReal BIon_h2o_h2o_offset;
      const static DoubleReal BIon_nh3_offset;
      const static DoubleReal BIon_nh3_h2o_offset;

      const static DoubleReal YIon_offset;
      const static DoubleReal YIon_h2o_offset;
      const static DoubleReal YIon_h2o_h2o_offset;
      const static DoubleReal YIon_nh3_offset;
      const static DoubleReal YIon_nh3_h2o_offset;

      const static DoubleReal AIon_offset;
      const static DoubleReal AIon_h2o_offset;
      const static DoubleReal AIon_nh3_offset;

      const static DoubleReal CIon_offset;

      const static DoubleReal ZIon_offset;

      ///generate vector containing all possible k-mer combinations of a set of size n
      void static nChooseKCombinations(int n, int k, std::vector<std::vector<int> > &combinations);

      //convenient output method for ion_type
      String static type_to_str(ion_type t)
      {
        switch(t)
        {
          case BIon: {return "BIon"; break;}
          case BIon2: {return "BIon2"; break;}
          case BIon_h2o: {return "BIon_h2o"; break;}
          case BIon_nh3: {return "BIon_nh3"; break;}
          case BIon_h2o_h2o: {return "BIon_h2o_h2o"; break;}
          case BIon_nh3_h2o: {return "BIon_nh3_h2o"; break;}
          case YIon: {return "YIon"; break;}
          case YIon2: {return "YIon2"; break;}
          case YIon_h2o: {return "YIon_h2o"; break;}
          case YIon_nh3: {return "YIon_nh3"; break;}
          case YIon_h2o_h2o: {return "YIon_h2o_h2o"; break;}
          case YIon_nh3_h2o: {return "YIon_nh3_h2o"; break;}
          case AIon: {return "AIon"; break;}
          case AIon_nh3: {return "AIon_nh3"; break;}
          case AIon_h2o: {return "AIon_h2o"; break;}
          case CIon: {return "CIon"; break;}
          case ZIon: {return "ZIon"; break;}
          default: {return "invalid"; break;}
        }
      }

      //convenient output method for ion_type
      ion_type static str_to_type(String s)
      {
        if(s=="BIon") return BIon;
        else if(s=="BIon2") return BIon2;
        else if(s=="BIon_h2o") return BIon_h2o;
        else if(s=="BIon_nh3") return BIon_nh3;
        else if(s=="BIon_h2o_h2o") return BIon_h2o_h2o;
        else if(s=="BIon_nh3_h2o") return BIon_nh3_h2o;
        else if(s=="YIon") return YIon;
        else if(s=="YIon2") return YIon2;
        else if(s=="YIon_h2o") return YIon_h2o;
        else if(s=="YIon_nh3") return YIon_nh3;
        else if(s=="YIon_h2o_h2o") return YIon_h2o_h2o;
        else if(s=="YIon_nh3_h2o") return YIon_nh3_h2o;
        else if(s=="AIon") return AIon;
        else if(s=="AIon_nh3") return AIon_nh3;
        else if(s=="AIon_h2o") return AIon_h2o;
        else if(s=="CIon") return CIon;
        else if(s=="ZIon") return ZIon;
        else return Invalid;
      }

      //returns fragment mass for ion_typ t, given the peptide mass (prefix for N-terminal, suffix for C-terminal)
      DoubleReal static getFragmentMass(ion_type t, DoubleReal m, DoubleReal parent_mass)
      {
        parent_mass=parent_mass-Residue::getInternalToFullMonoWeight();
        switch(t)
        {
          case BIon: return m+BIon_offset;break;
          case BIon2: return (m+BIon_offset+proton_m)/2.0;break;
          case BIon_h2o: return m + BIon_h2o_offset;break;
          case BIon_nh3: return m + BIon_nh3_offset;break;
          case BIon_h2o_h2o: return m + BIon_h2o_h2o_offset;break;
          case BIon_nh3_h2o: return m + BIon_nh3_h2o_offset;break;
          case YIon: return parent_mass - m + YIon_offset;break;
          case YIon2: return (parent_mass - m + YIon_offset + proton_m)/2;break;
          case YIon_h2o: return parent_mass - m + YIon_h2o_offset;break;
          case YIon_nh3: return parent_mass - m + YIon_nh3_offset;break;
          case YIon_h2o_h2o: return parent_mass - m + YIon_h2o_h2o_offset;break;
          case YIon_nh3_h2o: return parent_mass - m + YIon_nh3_h2o_offset;break;
          case AIon: return m + AIon_offset;break;
          case AIon_nh3: return m + AIon_nh3_offset;break;
          case AIon_h2o: return m + AIon_h2o_offset;break;
          case CIon: return m + CIon_offset;break;
          case ZIon: return parent_mass - m + ZIon_offset;break;
          default: std::cout<<"error in getRealOffsetMass"<<std::endl; return -1; break;
        }
      }


      //returns fragment mass for ion_typ t, given the peptide mass (prefix for N-terminal, suffix for C-terminal)
      DoubleReal static getPrefixMass(ion_type t, DoubleReal m, DoubleReal parent_mass)
      {
        parent_mass=parent_mass-Residue::getInternalToFullMonoWeight();
        switch(t)
        {
          case BIon: return m-BIon_offset;break;
          case BIon2: return (2*m)-BIon_offset-proton_m;break;
          case BIon_h2o: return m - BIon_h2o_offset;break;
          case BIon_nh3: return m - BIon_nh3_offset;break;
          case BIon_h2o_h2o: return m - BIon_h2o_h2o_offset;break;
          case BIon_nh3_h2o: return m - BIon_nh3_h2o_offset;break;
          case YIon: return parent_mass - m + YIon_offset;break;
          case YIon2: return parent_mass - 2*m + proton_m + YIon_offset;break;
          case YIon_h2o: return parent_mass - m + YIon_h2o_offset;break;
          case YIon_nh3: return parent_mass - m + YIon_nh3_offset;break;
          case YIon_h2o_h2o: return parent_mass - m + YIon_h2o_h2o_offset;break;
          case YIon_nh3_h2o: return parent_mass - m + YIon_nh3_h2o_offset;break;
          case AIon: return m - AIon_offset;break;
          case AIon_nh3: return m - AIon_nh3_offset;break;
          case AIon_h2o: return m - AIon_h2o_offset;break;
          case CIon: return m - CIon_offset;break;
          case ZIon: return parent_mass - m + ZIon_offset;break;
          default: std::cout<<"error in getRealOffsetMass"<<std::endl; return -1; break;
        }
      }



    protected:
    	
      ///allowed mass shift
      DoubleReal delta_;
	
			///mass precision for the amino acid masses
			DoubleReal precision_;

      ///If true Leucin and Isoleucin are not discriminated
      bool join_isobarics_;

      ///If true Lys and Gln are not discriminated
      bool join_Q_K_;

      ///The maximum gap length = max number of AA per edge in spectrum graph
      UInt max_span_;      

      ///generate vector containing all possible k-mer combinations of a set of size n
     void static nChooseKCombinations(int n, int k, std::vector<std::vector<int> > &combinations, std::vector<int> &tmp_comb);



  }; //IdSetup



//------------------------------------------------------------------------------------------------
//--------------------------------------------DEFINITIONS-----------------------------------------
//------------------------------------------------------------------------------------------------

/**function template
  reading array A from binary file and give to input vector<bool> arr
  templated to be used for both vectors A and A_len
   */
  template<typename T>
  std::streampos IdSetup::readVectorFromBinFile(String file, std::streampos pos, std::vector<T> &arr, Int max_elem)
  {
    std::ifstream in;

    in.exceptions(std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit);

    try
    {
      in.open(file.c_str(), std::ios::binary);

      T tmp_T;

      size_t size;

      in.seekg(pos, std::ios_base::beg);
      in.read((char *) & size, sizeof (size_t));

      //std::cout<<"bools to read: "<<size<<std::endl;
      if (max_elem>-1)
      {
        size = std::min(size, (size_t) max_elem);
      }

      std::cout<<"size to reserve: "<<size<<std::endl;
      arr.reserve(size);

      for (UInt i = 0; i < size; ++i)
      {
        in.read((char *) & tmp_T, sizeof (tmp_T));
        std::cout<<"reading: "<<tmp_T<<std::endl;
        arr.push_back(tmp_T);
      }

      pos = in.tellg();

      in.close();

    }
    catch (std::ifstream::failure &e)
    {
      if (!in.is_open())
      {
        std::cerr << "error opening binary file: " << file << " for reading vector data" << std::endl;
        throw;
      }
      else
      {
        in.close();
        std::cerr << "error reading vector data from binary file: " << file << std::endl;
        throw;
      }
    }

    return pos;
  }



///** write array A of type vector<bool> into a binary file
//templated for use with both vectors A and A_len
//*/
template<typename T>
std::streampos IdSetup::writeVectorToBinFile(String file, std::streampos pos, const std::vector<T> &arr, Int max_elem)
{
  std::cout<<"Entere writeVectorToBinFile test"<<std::endl;

  std::ofstream out;

    out.exceptions(std::ofstream::eofbit | std::ofstream::failbit | std::ofstream::badbit );

    try
    {
      if(pos!=(std::streampos)0)
      {
        //append to file
        out.open(file.c_str(), std::ios::binary | std::ios::app);
      }
      else
      {
        //rewrite file
        out.open(file.c_str(), std::ios::binary | std::ios::trunc);
      }

      size_t size=arr.size();

      if(max_elem>-1)
      {
        size=std::min(size,(size_t) max_elem);
      }

      // first entry in file is as UInt giving the number of elements that follow
      out.write((char *) &size,sizeof(size_t));

        for(UInt i = 0; i < size; ++i)
        {
          T tmp=arr[i];
          out.write((const char *) &tmp,sizeof(tmp));
          std::cout<<"writing: "<<tmp<<std::endl;
        }

        pos=out.tellp();

        out.close();

    }
    catch(std::ofstream::failure &e)
    {
      if(!out.is_open())
      {
        std::cerr<<"error opening binary file: "<<file<<" for writing AAmap data"<<std::endl;
        throw;
      }
      else
      {
        out.close();
        std::cerr<<"error writing AAmap data into binary file: "<<file<<std::endl;
        throw;
      }
    }

    return pos;
}



}//namespace

#endif /* IDSETUP_H_ */
