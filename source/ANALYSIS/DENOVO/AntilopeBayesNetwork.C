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

#include <OpenMS/ANALYSIS/DENOVO/AntilopeBayesNetwork.h>

//using namespace std;

namespace OpenMS
{

  BayesianNetwork::BayesianNetwork():
      XMLHandler("","0.3"),
      XMLFile("/Users/andreott/Work/ANTILOPE/DENOVOID/data/XMLBIF_0_3.xsd", "0.3")
  {
  }

  BayesianNetwork::BayesianNetwork(const BayesianNetwork & rhs):
      XMLHandler("","0.3"),
      XMLFile(),
      nodes_(rhs.nodes_)
  {
  }

  BayesianNetwork & BayesianNetwork::operator=(const BayesianNetwork & rhs)
  {
    if(this != &rhs)
    {
      nodes_=rhs.nodes_;
    }
    return *this;
  }


  void BayesianNetwork::load(const String &filename)
  {
    nodes_.clear();
    in_var_body_=false;

    file_ = filename;
    parse_(filename,this);

    std::map<String, Node>::iterator n_it;
    for(n_it=nodes_.begin(); n_it!=nodes_.end(); ++n_it)
    {
      n_it->second.num_outcomes=n_it->second.outcomes.size();
      n_it->second.outcome=0;
    }
  }

	void BayesianNetwork::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = String(sm_.convert(qname)).trim();

		//START
		if (tag_ == "BIF")
		{
			//check file version against schema version
			String file_version="";

			optionalAttributeAsString_(file_version,attributes,"version");
			if (file_version=="") file_version="0.3"; //default version is 1.0
			if (file_version.toDouble()>version_.toDouble())
			{
				warning(LOAD, "The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
			}
		}
		else if(tag_=="VARIABLE")
		{
			in_var_body_=true;
		}
	}

  void BayesianNetwork::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ = String(sm_.convert(qname)).trim();

    if(tag_=="VARIABLE")
    {
      in_var_body_=false;
    }
    else if(tag_=="DEFINITION")
    {
    }
    tag_="";
    return;
  }

    void BayesianNetwork::characters(const XMLCh* const chars, const XMLSize_t  /*length*/)
    {
      String value = ((String)sm_.convert(chars)).trim();

      if (tag_ == "NAME" && in_var_body_)
      {
        tag_ = "";
        for_node_id_=value;
        nodes_[for_node_id_].identifier=value;
        return;
      }
      else if (tag_ == "NAME" && !in_var_body_)
      {
        name_=value;
        return;
      }
      else if(tag_ == "OUTCOME")
      {
        nodes_[for_node_id_].outcomes.push_back(value);
        return;
      }
      else if(tag_ == "FOR")
      {
        for_node_id_=value;
        return;
      }
      else if(tag_ == "GIVEN")
      {
        nodes_[for_node_id_].parents.push_back(value);
        return;
      }
      else if(tag_ == "TABLE")
      {
        DoubleReal prob;
        std::stringstream value_stream(value);
        while(value_stream>>prob)
        {
          nodes_[for_node_id_].probs.push_back(prob);
        }
        return;
      }
    }


    void BayesianNetwork::print()
    {
      std::cout<<"Network name: "<<name_<<std::endl;

      std::map<String, Node>::iterator n_it;
      for(n_it=nodes_.begin(); n_it!=nodes_.end();++n_it)
      {
        std::cout<<std::endl<<"------------------------------"<<std::endl;
        std::cout<<"Node "<<n_it->second.identifier<<std::endl;
        std::vector<String>::iterator p_it;
        std::cout<<std::endl<<"Parent nodes: ";
        for(p_it=n_it->second.parents.begin(); p_it!=n_it->second.parents.end(); ++p_it)
        {
          std::cout<<*p_it<<" : ";
        }

        std::cout<<std::endl<<"Outcomes: ";
        std::vector<String>::iterator o_it;
        for(o_it=n_it->second.outcomes.begin(); o_it!=n_it->second.outcomes.end(); ++o_it)
        {
          std::cout<<*o_it<<" : ";
        }


        std::cout<<std::endl<<"Probabilities: ";
        std::vector<DoubleReal>::iterator pr_it;
        for(pr_it=n_it->second.probs.begin(); pr_it!=n_it->second.probs.end(); ++pr_it)
        {
          std::cout<<*pr_it<<" : ";
        }
        std::cout<<std::endl;
      }
    }

    DoubleReal BayesianNetwork::nodeCondProb(String identifier)
    {
      //compute the index of where to look in the cond_prob table
      const Node & node = nodes_[identifier];

      Size index = node.outcome;
      Size factor = node.num_outcomes;

      std::vector<String>::const_reverse_iterator parent_it;
      for(parent_it = node.parents.rbegin(); parent_it != node.parents.rend(); ++parent_it)
      {
        const Node & tmp_node = nodes_[*parent_it];
        index += tmp_node.outcome * factor;
        factor *= tmp_node.num_outcomes;
      }
      //std::cout<<"INDEX: "<<index<<std::endl;
      return node.probs[index];
    }

    DoubleReal BayesianNetwork::getLikelihood()
    {
      DoubleReal likelihood = 1;
      std::map<String, Node>::const_iterator n_it;
      for(n_it=nodes_.begin(); n_it!=nodes_.end(); ++n_it)
      {
        likelihood *= nodeCondProb(n_it->second.identifier);
      }
      return likelihood;
    }

    void BayesianNetwork::setOutComeIndex(String identifier, Size index)
    {
      std::map<String,Node>::iterator it = nodes_.find(identifier);
      if(it != nodes_.end() && index<it->second.num_outcomes)
      {
        it->second.outcome = index;
      }
    }

    void BayesianNetwork::getIdentifiers(std::set<String>&identifiers)
    {
      identifiers.clear();
      std::map<String, Node>::const_iterator node_it;
      for(node_it=nodes_.begin(); node_it!=nodes_.end(); ++node_it)
      {
        identifiers.insert(node_it->first);
      }
    }

}//Namespace


