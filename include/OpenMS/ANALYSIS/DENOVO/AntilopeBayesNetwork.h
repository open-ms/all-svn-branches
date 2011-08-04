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

#ifndef BIFXMLFILE_H
#define BIFXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

namespace OpenMS
{
  class BayesianNetwork:
      protected Internal::XMLHandler,
  public Internal::XMLFile
  {
  protected:

    //---BEGIN NESTED CLASS Node---
    class Node
    {
    public:
      ///The identifier, name for each node
      String identifier;

      //The value of the node, one of the possible outcomes
      Size outcome;

      //The names of the possible states for the node
      std::vector<String> outcomes;

      //The number of possible outcomes (must equal the size of outcomes - is stored to allow for faster computations)
      Size num_outcomes;

      //The identifiers of the nodes parents
      std::vector<String> parents;

      //The conditional probability values for the node
      std::vector<DoubleReal>probs;
    };//---END OF NESTED CLASS Node---


  public:
    ///Constructor
    BayesianNetwork();

    ///Destructor
    virtual ~BayesianNetwork(){}

    ///Copy constructor
    BayesianNetwork(const BayesianNetwork &);

    ///assignment operator
    BayesianNetwork & operator=(const BayesianNetwork &);


    ///return the probabilty of the given outcomes for the network
    DoubleReal getTotalCondProb();

    ///set the node outcome values (evidence)
    void setOutComeIndex(String identifier, Size index);

    /**
      @brief Loads the Bayesian network from file into an instance of the BayesianNetwork class

      The information is read in and the information is stored in the
      corresponding variables

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename);

    ///return the network name
    String getName()
    {
      return name_;
    }

    //return the keys of the nodes
    void getIdentifiers(std::set<String>&identifiers);

    ///print out the network --Basically for debugging
    void print();

    ///return the likelihood for the given set of observations
    DoubleReal getLikelihood();    

    void updateMembers_()
    {
    }

  protected:
    ///map containing the nodes of the network with their identifier as key
    std::map<String, Node>nodes_;

    ///name of the network
    String name_;

    ///get the conditional probability of node
    DoubleReal nodeCondProb(String node_identifier);

    // ---XML Parsing functions and variables---

    // Docu in base class
    virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

    // Docu in base class
    virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

    //Docu in base class
    virtual void characters(const XMLCh* const chars, const XMLSize_t  /*length*/);

    bool in_var_body_;

    String tag_;

    String for_node_id_;
  };




}//namespace

#endif // BIFXMLFILE_H
