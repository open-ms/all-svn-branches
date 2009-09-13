#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <iostream>
#include <fstream>

typedef int32_t LONG;
typedef float sup;
typedef LONG TT;

using namespace OpenMS;
using namespace std;

template<class T> int loadBinArray(char *filename, vector<vector<T> >& data)
{
  ifstream is(filename, ios::in | ios::binary);  //open file stream
  if(!is)
  {
    // An error occurred!
    is.clear();
    cout << "bgwaaak!" << endl; // todo throw real error
  }

  UInt c; // read header information of array size (m x n)
  is.read(reinterpret_cast<char*>(&c), sizeof(UInt) );
  const UInt m = c;
  is.read(reinterpret_cast<char*>(&c), sizeof(UInt) );
  const UInt n = c;

  //~ int length;
  //~ streampos pos_type;
  //~ pos_type = is.tellg();
  //~ // get length of file:
  //~ is.seekg(0, ios::end);
  //~ length = is.tellg();
  //~ is.seekg(pos_type);
  //~ cout << length/(sizeof(T)) <<endl;

  T *buff = new T[n]; //read the file according to the header into a buffer
  vector< vector<T> > inread(m, vector<T>(n));
  for(Size rows = 0; rows < m; ++rows)
  {
    is.read(reinterpret_cast<char*>(buff), n*sizeof(T) );
    if(is.good())
    {
      for(Size cols = 0; cols < n; ++cols)
      {
        inread[rows][cols] = buff[cols];
      }
    }
    else
    {
      cerr << "Data is incorrect\n"; return EXIT_FAILURE;
    }
  }

  //~ c= 666;
  //~ is.read(reinterpret_cast<char*>(&c), sizeof(UInt) ); cerr << c << "§" << is.good() << endl;
  //~ is.read(reinterpret_cast<char*>(&c), sizeof(UInt) ); cerr << c << "§" << is.good() << endl << true << endl ;

  delete[] buff; //return all resources to where they belong
  is.close();
  data = inread;
  return EXIT_SUCCESS;
}

int loadPklBin(char *filename, PeakMap& data)
{
  ifstream is(filename, ios::in | ios::binary);  //open file stream
  if(!is)
  {
    // An error occurred!
    is.clear();
    cout << "bgwaaak!" << endl; // todo throw real error
  }

  UInt c; // read header information of map size
  is.read(reinterpret_cast<char*>(&c), sizeof(int32_t) );
  const UInt num_specs = c;
  int16_t *num_buff = new int16_t[num_specs]; //read the spectra header according to the file header into a buffer
  is.read(reinterpret_cast<char*>(num_buff), num_specs*sizeof(int16_t) );
  if(!is.good())
  {
    cerr << "Data is incorrect\n"; return EXIT_FAILURE;
  }
  for(Size i = 0; i < num_specs; ++i)
  {
    float *peak_buff = new float[ (2*num_buff[i]+2) ]; //read the spectra according to the file and spectrum header into a buffer
    is.read(reinterpret_cast<char*>(peak_buff), (2*num_buff[i]+2)*sizeof(float) );
    if(is.good())
    {
      PeakSpectrum tmp_spec;
      Precursor tmp_prec;
      tmp_prec.setIsolationWindowLowerOffset(peak_buff[0]);
      tmp_prec.setIsolationWindowUpperOffset(peak_buff[0]);
      if(peak_buff[1]<4)
      {
        tmp_prec.setCharge(peak_buff[1]);
      }
      //~ else
      //~ {
        //~ tmp_prec.setCharge(2);
      //~ }
      tmp_spec.getPrecursors().push_back(tmp_prec);

      for(Size j = 2; j < (2*num_buff[i]+2); ++j)
      {
        Peak1D tmp_peak;
        tmp_peak.setMZ(peak_buff[j]);
        ++j;
        tmp_peak.setIntensity(peak_buff[j]);
        tmp_spec.push_back(tmp_peak);
      }
      data.push_back(tmp_spec);
    }
    else
    {
      cerr << "Data is incorrect\n"; return EXIT_FAILURE;
    }
    delete[] peak_buff; //return all resources to where they belong
  }

  delete[] num_buff; //return all resources to where they belong
  is.close();
  return EXIT_SUCCESS;
}
/* function specSet = sn_load_pklbin(filename, convertToDouble)
% function specSet = sn_load_pklbin(filename, convertToDouble)

if nargin<2 convertToDouble=1; end;
fid=fopen(filename,'r');  if fid<=0 fprintf(1,'Error opening file %s!\n',filename); specSet=[]; return; end;

numSpecs = fread(fid,1,'int32');   specSet = cell(numSpecs,5);
numPeaks = fread(fid,numSpecs,'int16=>double');
for i=1:numSpecs
    if convertToDouble data = fread(fid,2*numPeaks(i)+2,'float32=>double');
    else data = fread(fid,2*numPeaks(i)+2,'float32'); end;
    data2 = reshape(data,2,numPeaks(i)+1)';
    specSet{i,2} = data2(2:numPeaks(i)+1,:);
    specSet{i,3} = data2(1,1);
    specSet{i,4} = round(double(data2(1,2)));
    specSet{i,5} = '';
end

fclose(fid); */

void getNetAnnotated(PeakMap& map, vector< vector< pair<Size,DoubleReal> > >& blists, vector<DoubleReal>& ei, vector<DoubleReal>& tp, DoubleReal& peak_tol, DoubleReal& min_ei, DoubleReal& min_tp, bool include_iso = true)
{
  for(Size map_index = 0; map_index < map.size(); ++map_index)
  {

    if(map[map_index].size()!=0 && map[map_index].getPeptideIdentifications().size()!=0)
    {
      if(map[map_index].getPeptideIdentifications().front().getHits().front().getSequence().size()==1) //hack: just get the first available id
      {
        continue;
      }

      //sum total intensity
      Real total_intens(0);
      for(PeakSpectrum::const_iterator cit = map[map_index].begin(); cit != map[map_index].end(); ++cit)
      {
        if(cit->getMZ() > 5)
        {
          total_intens += cit->getIntensity();
        }
      }
      total_intens == 0 ? total_intens = 1 : total_intens;

      //todo peaks mowing? eigentlich schon bei input
      ParentPeakMower filter;
      filter.filterPeakSpectrum(map[map_index]);

      AASequence id = map[map_index].getPeptideIdentifications().front().getHits().front().getSequence(); //hack: just get the first available id
      //~ id.setNTerminalModification("MOD:09999"); // todo set fixed mods

      //OPENMS-way
      //get the theoretical spectra
      PeakSpectrum bs,ys/* ,bisos,yisos */;
      RichPeakSpectrum rps;
      TheoreticalSpectrumGenerator generator;
      generator.addPeaks(rps, id, Residue::BIon, 1);
      for(RichPeakSpectrum::Iterator it = rps.begin(); it != rps.end(); ++it)
      {
        bs.push_back(static_cast<Peak1D>(*it));
      }
      rps.clear();
      generator.addPeaks(rps, id, Residue::YIon, 1);
      for(RichPeakSpectrum::Iterator it = rps.begin(); it != rps.end(); ++it)
      {
        ys.push_back(static_cast<Peak1D>(*it));
      }

      //get the matches
      SpectrumAlignment sa;
      vector< pair< Size,Size > > alignment_B,alignment_Y;
      sa.getSpectrumAlignment(alignment_B, bs, map[map_index]);

      sa.getSpectrumAlignment(alignment_Y, ys, map[map_index]);

      //blists
      vector< pair<Size, DoubleReal> > matches_and_scores;
      for(Size s = 0; s < alignment_B.size(); ++s)
      {
        matches_and_scores.push_back(pair<Size, DoubleReal>(alignment_B[s].first,map[map_index][alignment_B[s].second].getIntensity()/total_intens));
      }
      blists[map_index] = matches_and_scores;

      //ipcountsB
      vector<Size> matched_indices_Bdiff(alignment_B.size()-1);
      for(Size s = 0; s < matched_indices_Bdiff.size(); ++s)
      {
        matched_indices_Bdiff[s] = alignment_B[s+1].first - alignment_B[s].first;
      }

/*    todo ip_counts_B/Y no use for them here
      vector<Size> ip_count_B(6,0); //n-th cell contains #(distance of consecutive matches == n)
      for(Size s = 0; s < matched_indices_Bdiff.size(); ++s)
      {
        switch(matched_indices_Bdiff[s])
        {
          case 0: cout << "err" << endl;  break; // todo - throw a real error
          case 1: ++ip_count_B[0]; break;
          case 2: ++ip_count_B[1]; break;
          case 3: ++ip_count_B[2]; break;
          case 4: ++ip_count_B[3]; break;
          case 5: ++ip_count_B[4]; break;
          default: ++ip_count_B[5]; break;
        }
      }

      //analogue ipcountY
      vector<Size> matched_indices_Ydiff(alignment_Y.size()-1);
      for(Size s = 0; s < matched_indices_Ydiff.size(); ++s)
      {
        matched_indices_Ydiff[s] = alignment_Y[s+1].first - alignment_Y[s].first;
      }

      vector<Size> ip_count_Y(6,0); //n-th cell contains #(distance of consecutive matches == n)
      for(Size s = 0; s < matched_indices_Ydiff.size(); ++s)
      {
        switch(matched_indices_Ydiff[s])
        {
          case 0: cout << "err" << endl;  break; // todo - throw a real error
          case 1: ++ip_count_Y[0]; break;
          case 2: ++ip_count_Y[1]; break;
          case 3: ++ip_count_Y[2]; break;
          case 4: ++ip_count_Y[3]; break;
          case 5: ++ip_count_Y[4]; break;
          default: ++ip_count_Y[5]; break;
        }
      } */

      // todo percs!!
      DoubleReal percs[8];
      percs[1] = (DoubleReal)alignment_B.size()/(DoubleReal)bs.size();
      percs[3] = (DoubleReal)alignment_Y.size()/(DoubleReal)ys.size();
      percs[5] = 0.0;
      for(Size s = 0; s < alignment_B.size(); ++s)
      {
        percs[5] += map[map_index][alignment_B[s].second].getIntensity();
      }
      percs[5] /= total_intens;
      percs[7] = 0.0;
      for(Size s = 0; s < alignment_B.size(); ++s)
      {
        percs[7] += map[map_index][alignment_Y[s].second].getIntensity();
      }
      percs[7] /= total_intens;

      if(percs[7]>percs[5]) //todo reverse spectrum todo find out why
      {
        //todo recalc percs recalc blists
        //maybe not replace but add reversed spectrum? see pathproj!
      }

      //remove unreliable ids
      if((percs[5]+percs[7])<min_ei || max(percs[1],percs[3])<min_tp)
      {
        //remove id
        map[map_index].getPeptideIdentifications().clear();
      }

      //replace id'ed spectrum with theoretical spec. (intens. = 0 or blist-entry(scores))
      PeakSpectrum tmp_spec;
      for(PeakSpectrum::Iterator it = bs.begin(); it != bs.end(); ++it)
      {
        Peak1D tmp_peak = *it;
        //~ tmp_peak.getMZ()== blist_entry?:tmp_peak.setIntensity(blists[map_index]):tmp_peak.setIntensity(0);
        tmp_spec.push_back(tmp_peak);
      }
      map[map_index] = tmp_spec;

/*    copious matlab implementation:
      //calc massesB/Y/Biso/Biso
      vector<DoubleReal> masses;
      for(Size res_index  = 0; res_index < id.size(); ++res_index)
      {
        masses.push_back(id.getResidue(res_index).getMonoWeight(Residue::Internal)); //getWeight see table like ionsource.com is Residue::Internal - todo more accurate might be something else
      }
      vector<DoubleReal> masses_B, masses_Y, masses_B_iso, masses_Y_iso;
      DoubleReal cum_mass(0);
      for(vector<DoubleReal>::const_iterator mass_it = masses.begin(); mass_it != masses.end()-1; ++mass_it) //only the n-1 first! masses --> only the ions
      {
        cum_mass += *mass_it;
        masses_B.push_back(cum_mass +b_offset ); //matlab implementation adds element-wise a 0 offset -maybe to correct internal-weight of AA to b-ion weight
        masses_B_iso.push_back(cum_mass +b_profile-b_profile ); //matlab implementation adds and subtracts element-wise a profile
      }
      cum_mass=0;
      for(vector<DoubleReal>::const_reverse_iterator mass_rit  = masses.rbegin(); mass_rit != masses.rend() -1 ; ++mass_rit) //only the n-1 last! masses (inverse) --> only the ions
      {
        cum_mass += *mass_rit;
        masses_Y.push_back(cum_mass +y_offset ); //matlab implementation adds element-wise a 18,0106 offset -maybe to correct internal-weight of AA to b-ion weight
        masses_Y_iso.push_back(cum_mass +y_profile-y_profile ); //matlab implementation adds and subtracts element-wise a profile
      }

      //calc diffs
      vector< vector<DoubleReal> > diffs_B(masses_B.size(),vector<DoubleReal>(map[map_index].size(),-1));
      vector< vector<DoubleReal> > diffs_Biso(masses_B.size(),vector<DoubleReal>(map[map_index].size(),-1));
      for(Size n = 0; n <= masses_B.size();++n)
      {
        for(Size m = 0; m <= map[map_index].size();++m)
        {
          diffs_B[n][m] = abs(masses_B[n]-map[map_index][m].getMZ());
          diffs_Biso[n][m] = abs(masses_B_iso[n]-map[map_index][m].getMZ());
        }
      }

      vector<bool> matches_B(masses_B.size(),false);
      vector<bool> matches_Biso(masses_B.size(),false);
      for(Size n = 0; n <= diffs_B.size();++n)
      {
        (*(min_element(diffs_B[n].begin(),diffs_B[n].end())) <= peak_tol) ? matches_B[n] = true : matches_B[n];
        (*(min_element(diffs_Biso[n].begin(),diffs_Biso[n].end())) <= peak_tol) ? matches_Biso[n] = true : matches_Biso[n];
      }

      vector<Size> matched_indices_B,matched_indices_Biso;
      for(Size s = 0; s < matches_B.size(); ++s)
      {
        if(matches_B[s] == true)
        {
          matched_indices_B.push_back(s);
        }
        if(matches_Biso[s] == true && include_iso)
        {
          matched_indices_B.push_back(s); //att! from 0 on!
        }
      }

      vector<Size> matched_indices_Bdiff(matched_indices_B.size()-1);
      for(Size s = 0; s < matched_indices_Bdiff.size(); ++s)
      {
        matched_indices_Bdiff[s] = matched_indices_B[s+1] - matched_indices_B[s];
      }

      vector<Size> ip_count_B(6,0); //n-th cell contains #(distance of consecutive matches == n)
      for(Size s = 0; s < matched_indices_Bdiff.size(); ++s)
      {
        switch(matched_indices_Bdiff[s])
        {
          case 0: cout << "err" << endl;  break; // todo - throw a real error
          case 1: ++ip_count_B[0]; break;
          case 2: ++ip_count_B[1]; break;
          case 3: ++ip_count_B[2]; break;
          case 4: ++ip_count_B[3]; break;
          case 5: ++ip_count_B[4]; break;
          default: ++ip_count_B[5]; break;
        }
      }

      //from which diffs_? cols i need the resp. max val? from matched_indices_? they come! this corresp. to the reduced blists{s} = blists{s}(find(blists{s}(:,1)>0),:); first col
      //the resp. max val? /total intensity are the second col to blists{s}(:,2)= blists{s}(:,2)/totalIntensity;

      */
    }
  }
}


void starBuilder(PeakMap& pairs, vector< pair<Size,Size> >& pair_indices, vector< pair<Size,DoubleReal> >& modification_positions, DoubleReal peak_tol, DoubleReal pm_tol, DoubleReal epsilon, DoubleReal noise_prob)
{

}
Int main()
{

  //input 1
  char* name = "stars_indices.bin";
  vector< vector<TT> > inread;
  int res = loadBinArray<TT>(name, inread);
  if(!res)
  {
    //test output
    cout << "inread.size(): " << inread.size() << endl;
    for(Size rows = 0; rows < inread.size(); ++rows)
    {
      cout << "inread["<< rows <<"].size(): " << inread[rows].size() << endl;
      for(Size cols = 0; cols < inread[rows].size(); ++cols)
      {
         //~ cout << inread[rows][cols] << "-";
      }
    }
    cout << endl;
  }

  //input 2
  MzMLFile mzml;
  PeakMap map; //MSExperiment<Peak1D> see StandardTypes
  // convert MzXML to MzData
  mzml.load("spectra/test.mzML",map);

  //input 3
  TextFile pep_annots("spectra/peptides_tagsearch.txt", true); // todo try-catch
  for(TextFile::iterator it = pep_annots.begin()+1; it != pep_annots.end(); ++it)
  {
    vector<String> tmp;
    it->split('\t',tmp);

    Size index, protein_index;
    DoubleReal peptide_tag_score,	perc_peptide_score;

    if((tmp.size() == 6))
    {
      index = tmp[0].toInt()-1;// todo try-catch
      peptide_tag_score = tmp[2].toDouble();// todo try-catch
      perc_peptide_score = tmp[3].toDouble();// todo try-catch
      protein_index = tmp[4].toInt();// todo try-catch
      Size true_index = inread[0][index]-1;
      if(true_index < map.size())
      {
        PeptideHit actual_hit(perc_peptide_score,0,-1,tmp[1].trim());
        actual_hit.setProteinAccessions(vector<String>(1,tmp[5]));
        PeptideIdentification pid;
        pid.setScoreType("perc_peptide_score");
        pid.insertHit(actual_hit);
        map[true_index].getPeptideIdentifications().push_back(pid);
      }
      else
      {
        cout << "no good file " << endl; break;
      }
    }
    else
    {
      cout << "no good file " << endl; break;
    }
  }

  //test output
  //~ for(PeakMap::iterator it = map.begin(); it != map.end(); ++it)
  //~ {
    //~ for(vector<PeptideIdentification>::const_iterator it_pid = it->getPeptideIdentifications().begin(); it_pid != it->getPeptideIdentifications().end(); ++it_pid)
    //~ {
      //~ for(vector<PeptideHit>::const_iterator it_hit = it_pid->getHits().begin(); it_hit != it_pid->getHits().end(); ++it_hit)
      //~ {
        //~ cout << it_hit->getSequence() << " ";
      //~ }
      //~ cout << " - ";
    //~ }
    //~ cout << " \'\' " << endl;
  //~ }

  //~ vector<PeakSpectrum> blists(map.size());
  //~ vector<Real> ei(map.size()), tp(map.size());

  AASequence test = map[0].getPeptideIdentifications().front().getHits().front().getSequence();

  cout << test << endl;

  cout << test[Size(0)].getMonoWeight(Residue::Internal) << " & " << test[Size(0)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(0)].getMonoWeight(Residue::BIon) << " & " << test[Size(0)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(1)].getMonoWeight(Residue::Internal) << " & " << test[Size(1)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(1)].getMonoWeight(Residue::BIon) << " & " << test[Size(1)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(2)].getMonoWeight(Residue::Internal) << " & " << test[Size(2)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(2)].getMonoWeight(Residue::BIon) << " & " << test[Size(2)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(3)].getMonoWeight(Residue::Internal) << " & " << test[Size(3)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(3)].getMonoWeight(Residue::BIon) << " & " << test[Size(3)].getMonoWeight(Residue::YIon) << endl;

  test.setNTerminalModification("MOD:09999");

    cout << test[Size(0)].getMonoWeight(Residue::Internal) << " & " << test[Size(0)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(0)].getMonoWeight(Residue::BIon) << " & " << test[Size(0)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(1)].getMonoWeight(Residue::Internal) << " & " << test[Size(1)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(1)].getMonoWeight(Residue::BIon) << " & " << test[Size(1)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(2)].getMonoWeight(Residue::Internal) << " & " << test[Size(2)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(2)].getMonoWeight(Residue::BIon) << " & " << test[Size(2)].getMonoWeight(Residue::YIon) << endl;
  cout << test[Size(3)].getMonoWeight(Residue::Internal) << " & " << test[Size(3)].getMonoWeight(Residue::NTerminal ) << " & " << test[Size(3)].getMonoWeight(Residue::BIon) << " & " << test[Size(3)].getMonoWeight(Residue::YIon) << endl;


      vector<bool> matches_B(10,false);
      matches_B[0]=true;matches_B[4]=true;matches_B[9]=true;
      vector<Size> matched_indices;
      for(Size s = 0; s < matches_B.size(); ++s)
      {
        if(matches_B[s] == true){ matched_indices.push_back(s); }
      }
      for(Size s = 0; s < matches_B.size(); ++s)
      {
        cout << matches_B[s];
      }
      cout << endl;

      vector<Size> matched_indices_(matched_indices.size()-1);
      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        matched_indices_[s] = matched_indices[s+1] - matched_indices[s];
      }

      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        cout << matched_indices_[s];
      }
      cout << endl;

      vector<Size> ip_count_B(6,0); //n-te zelle enthält #(abstand consekutiver matches == n)
      for(Size s = 0; s < matched_indices_.size(); ++s)
      {
        cout << matched_indices_[s] << endl;
        switch(matched_indices_[s])
        {
          case 0: cout << "err" << endl; break;
          case 1: ++ip_count_B[0]; break;
          case 2: ++ip_count_B[1]; break;
          case 3: ++ip_count_B[2]; break;
          case 4: ++ip_count_B[3]; break;
          case 5: ++ip_count_B[4]; break;
          default: ++ip_count_B[5]; break;
        }
      }

      for(Size s = 0; s < ip_count_B.size(); ++s)
      {
        cout << ip_count_B[s];
      }
      cout << endl;

      //get the theoretical spectra
      test = "YYTS";
      PeakSpectrum bs,ys,bisos,yisos;
      RichPeakSpectrum rps;
      TheoreticalSpectrumGenerator generator;
      generator.addPeaks(rps, test, Residue::BIon, 1);

      cout << test << endl;
      for(Size i = 0; i < rps.size(); ++i)
      {
        cout << static_cast<Peak1D>(rps[i]) << endl;
      }

      SpectrumAlignment sa;
      vector< pair< Size,Size > > alignment;
      //~ sa.getSpectrumAlignment(alignment, rps, map[0]);


      rps.clear();
      generator.addPeaks(rps, test, Residue::YIon, 1);

      cout << test << endl;
      for(Size i = 0; i < rps.size(); ++i)
      {
        cout << static_cast<Peak1D>(rps[i]) << endl;
      }

      test = "S";
      cout << test.getMonoWeight(Residue::YIon,0) << endl;
      cout << test.getMonoWeight(Residue::CTerminal,0) << endl;
      cout << test.getMonoWeight(Residue::YIon,1) << endl;
      cout << test.getMonoWeight(Residue::CTerminal,1) << endl;
      cout << test.getMonoWeight(Residue::YIon,2) << endl;
      cout << test.getMonoWeight(Residue::CTerminal,2) << endl;
      cout << test.getMonoWeight(Residue::Internal,2) << endl;
      test = "Y";
      cout << test.getMonoWeight(Residue::BIon,0) << endl;
      cout << test.getMonoWeight(Residue::NTerminal,0) << endl;
      cout << test.getMonoWeight(Residue::BIon,1) << endl;
      cout << test.getMonoWeight(Residue::NTerminal,1) << endl;
      cout << test.getMonoWeight(Residue::BIon,2) << endl;
      cout << test.getMonoWeight(Residue::NTerminal,2) << endl;
      cout << test.getMonoWeight(Residue::Internal,2) << endl;


  //read a binarray int32_t
  char* name_1 = "kept_indices.bin";
  vector< vector<int32_t> > inread_1;
  int res_ = loadBinArray<int32_t>(name_1, inread_1);
  if(!res_)
  {
    cout << name_1 << endl;
  }

  //input test
  char* name_2 = "aligns/pairs_modpos.bin";
  vector< vector<sup> > inread_2;
  res_ = loadBinArray<sup>(name_2, inread_2);
  if(!res_)
  {
    cout << name_2 << endl;
    //test output
    cout << "inread.size(): " << inread_2.size() << endl;
    for(Size rows = 0; rows < inread_2.size(); ++rows)
    {
      cout << "inread["<< rows <<"].size(): " << inread_2[rows].size() << endl;
      for(Size cols = 0; cols < inread_2[rows].size(); ++cols)
      {
         cout << inread[rows][cols] << "-";
      }
    }
    cout << endl;
  }

  //~ //load a mgf
  //~ String filename_pairs("pairs_pklbin.mgf");
  PeakMap pairs_map;
  //~ FileHandler fh;
  //~ FileTypes::Type file_type = fh.getType(filename_pairs);
  //~ if (file_type!=FileTypes::MGF)
  //~ {
    //~ cout << "Open file error: Could not determine file type of '"<<filename_pairs<<"'!";
    //~ return EXIT_FAILURE;
  //~ }
  //~ fh.loadExperiment(filename_pairs, pairs_map, file_type);
  char* pairs_file = "aligns/pairs.pklbin";
  loadPklBin(pairs_file, pairs_map);



  return EXIT_SUCCESS;
} //end of main
