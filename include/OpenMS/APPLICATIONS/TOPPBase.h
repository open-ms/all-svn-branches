// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors:  Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_TOPPBASE_H
#define OPENMS_APPLICATIONS_TOPPBASE_H

#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/IntList.h>
#include <OpenMS/DATASTRUCTURES/DoubleList.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/DocumentIDTagger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>


#include <iostream>
#include <fstream>
#include <limits>

class QStringList;


namespace OpenMS
{
	class ConsensusMap;

  namespace Exception
  {
    /// An unregistered parameter was accessed
    class OPENMS_DLLAPI UnregisteredParameter
          : public Exception::BaseException
    {
      public:
        UnregisteredParameter( const char* file, int line, const char* function, const String& parameter )
            : BaseException( file, line, function, "UnregisteredParameter", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
    /// A parameter was accessed with the wrong type
    class OPENMS_DLLAPI WrongParameterType
          : public Exception::BaseException
    {
      public:
        WrongParameterType( const char* file, int line, const char* function, const String& parameter )
            : BaseException( file, line, function, "WrongParameterType", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
    /// A required parameter was not given
    class OPENMS_DLLAPI RequiredParameterNotGiven
          : public Exception::BaseException
    {
      public:
        RequiredParameterNotGiven( const char* file, int line, const char* function, const String& parameter )
            : BaseException( file, line, function, "RequiredParameterNotGiven", parameter )
        {
          globalHandler.setMessage( what_ );
        }
    };
  }

  /**
  	 @brief Base class for TOPP applications.

  	This base class implements functionality used in most TOPP tools:
  	- parameter handling
  	- file handling
  	- progress logging

  	If you want to create a new TOPP tool, please take care of the following:
  	- derive a new class from this class
  	- implement the registerOptionsAndFlags_ and main_ methods
  	- add a doxygen page for the tool and add the page to TOPP.doxygen
  	- hide the derived class in the OpenMS documentation by using doxygen condition macros.
  
		@todo write subsections if type was given or if no type is used; see MascotAdapterOnline (Andreas)

    @todo: replace writeLog_, writeDebug_ with a logger concept
           we'd need something like -VLevels <LOGGERS> to specify which loggers shall print something
           the '-log' flag should clone all output to the log-file (maybe with custom <LOGGERS), which can either be specified directly or is
              equal to '-out' (if present) with a ".log" suffix
           maybe a new LOGGER type (TOPP), which is only useable on TOPP level?


  */
  class OPENMS_DLLAPI TOPPBase
  {
    public:

      /// Exit codes
      enum ExitCodes
      {
        EXECUTION_OK,
        INPUT_FILE_NOT_FOUND,
        INPUT_FILE_NOT_READABLE,
        INPUT_FILE_CORRUPT,
        INPUT_FILE_EMPTY,
        CANNOT_WRITE_OUTPUT_FILE,
        ILLEGAL_PARAMETERS,
        MISSING_PARAMETERS,
        UNKNOWN_ERROR,
        EXTERNAL_PROGRAM_ERROR,
        PARSE_ERROR,
        INCOMPATIBLE_INPUT_DATA,
        INTERNAL_ERROR,
        UNEXPECTED_RESULT
    	};

      /**
      	@brief Constructor
      	
      	@param name Tool name.
      	@param description Short description of the tool (one line).
      	@param official If this is an official TOPP tool contained in the OpenMS/TOPP release.
      	                If @em true the tool name is checked against the list of TOPP tools and a warning printed if missing.
				@param id_tag_support Does the TOPP tool support unique DocumentIdentifier assignment?! The default is false.
							 In the default case you cannot use the -id_pool argument when calling the TOPP tool (it will terminate during init)
      	@param version Optional version of the tools (if empty, the version of OpenMS/TOPP is used).
      */
      TOPPBase(const String& name, const String& description, bool official=true, bool id_tag_support=false, const String& version="");

      /// Destructor
      virtual ~TOPPBase();

			/// @todo to be documented (Chris, Andreas)
			void checkTOPPIniFile(const String& tool_path);

      /// Main routine of all TOPP applications
      ExitCodes main(int argc, const char** argv);
			
      ///Stuct that captures all information of a command line parameter
      struct ParameterInformation
      {
        /// Parameter types
        enum ParameterTypes
        {
          NONE = 0,       ///< Undefined type
          STRING,         ///< String parameter
          INPUT_FILE,			///< String parameter that denotes an input file
          OUTPUT_FILE,    ///< String parameter that denotes an output file
          DOUBLE,         ///< Floating point number parameter
          INT,            ///< Integer parameter
          STRINGLIST,     ///< More than one String Parameter
          INTLIST,        ///< More than one Integer Parameter
          DOUBLELIST,     ///< More than one String Parameter
          INPUT_FILE_LIST,///< More than one String Parameter that denotes input files
          OUTPUT_FILE_LIST,///< More than one String Parameter that denotes output files
          FLAG,           ///< Parameter without argument
          TEXT,           ///< Left aligned text, see addText_
          NEWLINE					///< An empty line, see addEmptyLine_
        };

        /// name of the parameter (internal and external)
        String name;
        /// type of the parameter
        ParameterTypes type;
        /// default value of the parameter stored as string
        DataValue default_value;
        /// description of the parameter
        String description;
        /// argument in the description
        String argument;
        /// flag that indicates if this parameter is required i.e. it must differ from the default value
        bool required;
        /// flag the indicates that the parameter is advanced (this is used for writing the INI file only)
        bool advanced;
        /// StringList for special tags
        StringList tags;

				///@name Restrictions for different parameter types
				//@{
				std::vector<String> valid_strings;
				Int min_int;
				Int max_int;
				DoubleReal min_float;
				DoubleReal max_float;
				//@}
				
        /// Constructor that takes all members in declaration order
        ParameterInformation( const String& n, ParameterTypes t, const String& arg, const DataValue& def, const String& desc, bool req, bool adv, const StringList& tag_values = StringList() )
        	: name(n),
	          type(t),
	          default_value(def),
	          description(desc),
	          argument(arg),
	          required(req),
	          advanced(adv),
            tags(tag_values),
	          valid_strings(),
            min_int(-std::numeric_limits<Int>::max()),
        	  max_int(std::numeric_limits<Int>::max()),
        	  min_float(-std::numeric_limits<DoubleReal>::max()),
        	  max_float(std::numeric_limits<DoubleReal>::max())
        {
        }

        ParameterInformation()
          : name(),
            type( NONE ),
            default_value(),
            description(),
            argument(),
            required(true),
        	  advanced(false),
            tags(),
            valid_strings(),
            min_int(-std::numeric_limits<Int>::max()),
        	  max_int(std::numeric_limits<Int>::max()),
        	  min_float(-std::numeric_limits<DoubleReal>::max()),
        	  max_float(std::numeric_limits<DoubleReal>::max())
        {
        }

        ParameterInformation& operator=( const ParameterInformation& rhs )
        {
          if ( &rhs == this ) return *this;

          name = rhs.name;
          type = rhs.type;
          default_value = rhs.default_value;
          description = rhs.description;
          argument = rhs.argument;
          required = rhs.required;
          advanced = rhs.advanced;
          tags == rhs.tags;
          valid_strings = rhs.valid_strings;
          min_int = rhs.min_int;
          max_int = rhs.max_int;
          min_float = rhs.min_float;
          max_float = rhs.max_float;
          
          return *this;
        }
      };


    private:

      /// Tool name.  This is assigned once and for all in the constructor.
      String const tool_name_;

      /// Tool description. This is assigned once and for all in the constructor.
      String const tool_description_;

			/// Tool indicates it supports assignment of unique DocumentID from IDPool
			bool id_tag_support_;

			/// Instance of DocumentIDTagger, which can be accessed using getDocumentIDTagger_()
			DocumentIDTagger id_tagger_;

      ///Instance number
      Int const instance_number_;

      ///Location in the ini file where to look for parameters.
      String const ini_location_;

      /// No default constructor.  It is "declared away".
      TOPPBase();

      /// No default copy constructor.  It is "declared away".
      TOPPBase( const TOPPBase& );

      /// All parameters relevant to this invocation of the program.
      Param param_;

      /// All parameters specified in the ini file
      Param param_inifile_;

      /// Parameters from command line
      Param param_cmdline_;

      /// Parameters from instance section
      Param param_instance_;

      /// Parameters from common section with tool name.
      Param param_common_tool_;

      /// Parameters from common section without tool name.
      Param param_common_;

      ///Log file stream.  Use the writeLog_() and writeDebug_() methods to access it.
      mutable std::ofstream log_;

      /** @brief Ensures that at least some default logging destination is
      		opened for writing in append mode.

      		@note This might be invoked at various places early in the startup
      		process of the TOPP tool.  Thus we cannot consider the ini file here.
      		The final logging destination is determined in main().
      */
      void enableLogging_() const;

      /// Storage location for parameter information
      std::vector<ParameterInformation> parameters_;

      /**
      	@brief This method should return the default parameters for subsections.

      	It is called once for each registered subsection, when writing the an example ini file.

      	Reimplement this method to set the defaults written in the 'write_ini' method.
      	
      	@note Make sure to set the 'advanced' flag of the parameters right in order to hide certain parameters from unexperienced users.
      */
      virtual Param getSubsectionDefaults_( const String& section ) const;

      /// Storage location and description for allowed subsections
      std::map<String, String> subsections_;

      /// Storage location and description for allowed subsections from TOPP tool's command-line parameters
      std::map<String, String> subsections_TOPP_;

      /**
      	@name Internal parameter handling
       */
      //@{
      /**
      	 @brief Return the value of parameter @p key as a string or @p default_value if this value is not set.
      
      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      String getParamAsString_( const String& key, const String& default_value = "" ) const;

      /**
      	 @brief Return the value of parameter @p key as an integer or @p default_value if this value is not set.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      Int getParamAsInt_( const String& key, Int default_value = 0 ) const;

      /**
      	 @brief Return the value of parameter @p key as a double or @p default_value if this value is not set.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      double getParamAsDouble_( const String& key, double default_value = 0 ) const;
      
      /**
         @brief Return the value of parameter @p key as a StringList or @p default_value if this value is not set
         
         @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      StringList getParamAsStringList_(const String& key,const StringList& default_value) const;
      
      /**
         @brief Return the value of parameter @p key as a IntList or @p default_value if this value is not set
         
         @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      IntList getParamAsIntList_(const String& key,const IntList& default_value) const;

      /**
         @brief Return the value of parameter @p key as a DoubleList or @p default_value if this value is not set
         
         @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      DoubleList getParamAsDoubleList_(const String& key,const DoubleList& default_value) const;
      
      /**
      	 @brief Return the value of flag parameter @p key as bool.

      	 Only the string values 'true' and 'false' are interpreted.
      	 
      	 @exception Exception::InvalidParameter is thrown for non-string parameters and string parameters with values other than 'true' and 'false'.

      	 @note See getParam_(const String&) const for the order in which parameters are searched.
      */
      bool getParamAsBool_( const String& key) const;
      
      /**
      	 @brief Return the value @p key of parameters as DataValue. DataValue::EMPTY indicates that a parameter was not found.

      	 Parameters are searched in this order:
      	 -# command line
      	 -# instance section, e.g. "TOPPTool:1:some_key", see getIniLocation_().
      	 -# common section with tool name,  e.g. "common:ToolName:some_key"
      	 -# common section without tool name,  e.g. "common:some_key"
      	 
      	 where "some_key" == key in the examples.
      */
      const DataValue& getParam_( const String& key ) const;
      
      /// Returns the default parameters
      Param getDefaultParameters_() const;

			/// Returns the user defaults for the given tool, if any default parameters are stored in the users home
			Param getToolUserDefaults_(const String& tool_name) const;
      //@}

    protected:
			///Version string (if empty, the OpenMS/TOPP version is printed)
			String version_;
			
      /**
      	@brief Returns the location of the ini file where parameters are taken
      	from.  E.g. if the command line was <code>TOPPTool -instance 17</code>, then
      	this will be <code>"TOPPTool:17:"</code>.  Note the ':' at the end.

      	This is assigned during tool startup, depending on the command line but (of course) not depending on ini files.
      */
      const String& getIniLocation_() const
      {
        return ini_location_;
      }
			
			///Returns the tool name
			const String& toolName_() const;
			
      /**
      	@name Parameter handling

      	Use the methods registerStringOption_, registerInputFile_, registerOutputFile_, registerDoubleOption_, 
      	registerIntOption_ and registerFlag_ in order to register parameters in registerOptionsAndFlags_.

      	To access the values of registered parameters in the main_ method use methods
      	getStringOption_ (also for input and output files), getDoubleOption_, getIntOption_,getStringList_(also for input and output file lists),getIntList_,getDoubleList_, and getFlag_.

				The values of certain options can be restricted using: setMinInt_, setMaxInt_, setMinFloat_, 
				setMaxFloat_, setValidStrings_ and setValidFormats_.

      	In order to format the help output, the methods addEmptyLine_ and addText_ can be used.
      */
      //@{
      /**
      	 @brief Sets the valid command line options (with argument) and flags (without argument).

      	 The options '-ini' '-log' '-instance' '-debug' and the flag '--help' are automatically registered.
      */
      virtual void registerOptionsAndFlags_() = 0;

      /**
      	@brief Registers a string option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerStringOption_(const String& name, const String& argument, const String& default_value, const String& description, bool required = true, bool advanced = false);
			
			/**
				@brief Sets the valid strings for a string option or a whole string list
				
				@exception Exception::ElementNotFound is thrown if the parameter is unset or not a string parameter
				@exception Exception::InvalidParameter is thrown if the valid strings contain comma characters
			*/
			void setValidStrings_(const String& name, const std::vector<String>& strings);
			
      /**
      	@brief Registers an input file option.
				
				Input files behave like string options, but are automatically checked with inputFileReadable_()
				when the option is accessed in the TOPP tool.
				
      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      	@param advanced A list of tags, e.g. 'skipexists', specifying the handling of the input file (e.g. when its an executable)
                        Valid tags: 'skipexists' - will prevent checking if the given file really exists (useful for an executable in global PATH)
      */
      void registerInputFile_( const String& name, const String& argument, const String& default_value, const String& description, bool required = true, bool advanced = false, const StringList& tags=StringList() );

      /**
      	@brief Registers an output file option.
				
				Output files behave like string options, but are automatically checked with outputFileWritable_()
				when the option is accessed in the TOPP tool.
				
      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerOutputFile_( const String& name, const String& argument, const String& default_value, const String& description, bool required = true, bool advanced = false );

			/**
				@brief Sets the formats for a input/output file option or for all members of an input/output file lists
				
				Setting the formats causes a check for the right file format (input file) or the right file extension (output file).
				This check is performed only, when the option is accessed in the TOPP tool.
				When @p force_OpenMS_format is set, only formats known to OpenMS internally are allowed (default).

				@exception Exception::ElementNotFound is thrown if the parameter is unset or not a file parameter
				@exception Exception::InvalidParameter is thrown if an unknown format name is used (@see FileHandler::Type)
			*/
			void setValidFormats_(const String& name, const std::vector<String>& formats, const bool force_OpenMS_format=true);
			

      /**
      	@brief Registers a double option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerDoubleOption_( const String& name, const String& argument, double default_value, const String& description, bool required = true, bool advanced = false );

			/**
				@brief Sets the minimum value for the integer parameter(can be a list of integers,too) @p name. 
				
				@exception Exception::ElementNotFound is thrown if @p name is not found or if the parameter type is wrong
			*/			
			void setMinInt_(const String& name, Int min);
			/**
				@brief Sets the maximum value for the integer parameter(can be a list of integers,too) @p name. 
				
					@exception Exception::ElementNotFound is thrown if @p name is not found or if the parameter type is wrong
		*/
			void setMaxInt_(const String& name, Int max);
			/**
				@brief Sets the minimum value for the floating point parameter(can be a list of floating points,too) @p name. 
				
				@exception Exception::ElementNotFound is thrown if @p name is not found or if the parameter type is wrong
			*/
			void setMinFloat_(const String& name, DoubleReal min);
			/**
				@brief Sets the maximum value for the floating point parameter(can be a list of floating points,too) @p name. 
				
				@exception Exception::ElementNotFound is thrown if @p name is not found or if the parameter type is wrong
			*/
			void setMaxFloat_(const String& name, DoubleReal max);

      /**
      	@brief Registers an integer option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerIntOption_(const String& name, const String& argument, 
															Int default_value, const String& description, 
															bool required = true, bool advanced = false );

      /**
      	@brief Registers a list of integers option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.

      */
      void registerIntList_( const String& name, const String& argument, IntList default_value, const String& description, bool required = true, bool advanced = false );
      
     /**
      	@brief Registers a list of doubles option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerDoubleList_( const String& name, const String& argument, DoubleList default_value, const String& description, bool required = true, bool advanced = false );
      
     /**
      	@brief Registers a list of strings option.

      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerStringList_( const String& name, const String& argument, StringList default_value, const String& description, bool required = true, bool advanced = false );

     /**
      	@brief Registers a list of input files option.
        
        A list of input files behaves like a StringList, but are automatically checked with inputFileWritable_()
				when the option is accessed in the TOPP tool.
        
      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerInputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required = true, bool advanced = false );

     /**
      	@brief Registers a list of output files option.
        
        A list of output files behaves like a StringList, but are automatically checked with outputFileWritable_()
				when the option is accessed in the TOPP tool.
        
      	@param name Name of the option in the command line and the INI file
      	@param argument Argument description text for the help output
      	@param default_value Default argument
      	@param description Description of the parameter. Indentation of newline is done automatically.
      	@param required If the user has to provide a value i.e. if the value has to differ from the default (checked in get-method)
      	@param advanced If @em true, this parameter is advanced and by default hidden in the GUI.
      */
      void registerOutputFileList_( const String& name, const String& argument, StringList default_value, const String& description, bool required = true, bool advanced = false );
      
      /// Registers a flag
      void registerFlag_( const String& name, const String& description, bool advanced = false );

      /**
      	@brief Registers an allowed subsection in the INI file (usually from OpenMS algorithms).

      	Use this method to register subsections that are passed to algorithms.

      	@see checkParam_
      */
      void registerSubsection_( const String& name, const String& description );
      
      /**
      	@brief Registers an allowed subsection in the INI file originating from the TOPP tool itself.

      	Use this method to register subsections which is created by a commandline param (registered by e.g. registerDoubleOption_() )
        and contains a ':' in its name. This is done to distinguish these parameters from normal subsections,
        which are filled by calling 'getSubsectionDefaults_()'. This is not necessary for here.

      	@see checkParam_
      */
      void registerTOPPSubsection_( const String& name, const String& description );


      /// Adds an empty line between registered variables in the documentation.
      void addEmptyLine_();

      /// Adds a left aligned text between registered variables in the documentation e.g. for subdividing the documentation.
      void addText_( const String& text );

      /**
        @brief Returns the value of a previously registered string option

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
        @exception Exception::RequiredParameterNotGiven is if a required parameter is not present
        @exception Exception::WrongParameterType is thrown if the parameter has the wrong type
        @exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
      */
      String getStringOption_( const String& name ) const;

      /**
      	@brief Returns the value of a previously registered double option

				@exception Exception::UnregisteredParameter is thrown if the parameter was not registered
				@exception Exception::RequiredParameterNotGiven is if a required parameter is not present
				@exception Exception::WrongParameterType is thrown if the parameter has the wrong type
				@exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
			*/
			DoubleReal getDoubleOption_( const String& name ) const;

      /**
      	@brief Returns the value of a previously registered integer option

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
        @exception Exception::RequiredParameterNotGiven is if a required parameter is not present
        @exception Exception::WrongParameterType is thrown if the parameter has the wrong type
        @exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
      */
      Int getIntOption_( const String& name ) const;
      
      /**
      	@brief Returns the value of a previously registered StringList

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
        @exception Exception::RequiredParameterNotGiven is if a required parameter is not present
        @exception Exception::WrongParameterType is thrown if the parameter has the wrong type
        @exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
      */
      StringList getStringList_( const String& name ) const;

      /**
      	@brief Returns the value of a previously registered IntList

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
        @exception Exception::RequiredParameterNotGiven is if a required parameter is not present
        @exception Exception::WrongParameterType is thrown if the parameter has the wrong type
        @exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
      */
      IntList getIntList_( const String& name ) const;
      
      /**
      	@brief Returns the value of a previously registered DoubleList

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
        @exception Exception::RequiredParameterNotGiven is if a required parameter is not present
        @exception Exception::WrongParameterType is thrown if the parameter has the wrong type
        @exception Exception::InvalidParameter is thrown if the parameter restrictions are not met
      */
      DoubleList getDoubleList_( const String& name ) const;
      
      ///Returns the value of a previously registered flag
      bool getFlag_( const String& name ) const;

      /**
        @brief Finds the entry in the parameters_ array that has the name @p name

        @exception Exception::UnregisteredParameter is thrown if the parameter was not registered
      */
      const ParameterInformation& findEntry_( const String& name ) const;

      /**
      	@brief Return <em>all</em> parameters relevant to this TOPP tool.

      	Returns a Param that contains everything you can get by the getParamAs...() methods.
      */
      Param const& getParam_() const;

      /**
      	@brief Checks top-level entries of @p param according to the information during registration

      	Only top-level entries and allowed subsections are checked.
      	Checking the content of the subsection is the duty of the algorithm it is passed to.

      	This method does not abort execution of the tool, but will warn the user through stderr!
      	It is called automatically in the main method.

      	@param param Parameters to check
      	@param filename The source file name
      	@param location Exact location inside the source file
      */
      void checkParam_( const Param& param, const String& filename, const String& location ) const;
      //@}

      /// make a string console friendly
      String breakString_(const String& input, const Size line_len, const Size indentation, const Size max_lines) const;
      
      /// read console settings for output shaping
      void readConsoleSize_();

      /// Prints the tool-specific command line options and appends the common options.
      void printUsage_();

      /// The actual "main" method.  main_() is invoked by main().
      virtual ExitCodes main_(int argc , const char** argv) = 0;

      ///@name Debug and Log output
      //@{
      /// Writes a string to the log file and to std::cout
      void writeLog_( const String& text ) const;

      /// Writes a @p text to the log file and to std::cout if the debug level is at least @p min_level
      void writeDebug_( const String& text, UInt min_level ) const;

      /// Writes a String followed by a Param to the log file and to std::cout if the debug level is at least @p min_level
      void writeDebug_( const String& text, const Param& param, UInt min_level ) const;
      //@}


      /**
      	@name File IO checking methods

      	Methods used to check the validity of input and output files in main_.
				
				Checking input and output files is only necessary, if you did register the file as string option,
				e.g. when only a file prefix is given which is completed in the program.	
				
      	The exceptions thrown in these methods are catched in the main method of this class.
      	They do not have to be handled in the tool itself!
      */
      //@{
      /**
        @brief Checks if an input file exists, is readable and is not empty

        The @em filename is a URI to the file to be read and @em param_name gives the name of the parameter
        , e.g. "in" which specified the filename (this is useful for error messages when the file cannot be read, so the
        user can immediately see which parameter to change). If no parameter is responsible for the
        name of the input file, then leave @em param_name empty.

        @exception Exception::FileNotFound is thrown if the file is not found
        @exception Exception::FileNotReadable is thrown if the file is not readable
        @exception Exception::FileEmpty is thrown if the file is empty
      */
      void inputFileReadable_( const String& filename, const String& param_name) const;

      /**
        @brief Checks if an output file is writeable

        The @em filename is a URI to the file to be written and @em param_name gives the name of the parameter
        , e.g. "out" which specified the filename (this is useful for error messages when the file cannot be written, so the
        user can immediately see which parameter to change). If no parameter is responsible for the
        name of the output file, then leave @em param_name empty.

        @exception Exception::UnableToCreateFile is thrown if the file cannot be created
      */
      void outputFileWritable_( const String& filename, const String& param_name) const;
      //@}

      /// Helper function that parses a range string ([a]:[b]) into two variables
      void parseRange_( const String& text, double& low, double& high ) const;

      ///Type of progress logging
      ProgressLogger::LogType log_type_;

			///@name Data processing auxilary functions
      //@{
      
      ///Data processing setter for consensus maps
      void addDataProcessing_(ConsensusMap& map, const DataProcessing& dp) const;

      ///Data processing setter for feature maps
      template<typename FeatureType>
      void addDataProcessing_(FeatureMap<FeatureType>& map, const DataProcessing& dp) const
      {
        map.getDataProcessing().push_back(dp);          
      }
			
			///Data processing setter for peak maps
      template<typename PeakType>
      void addDataProcessing_(MSExperiment<PeakType>& map, const DataProcessing& dp) const
      {
        for (Size i=0; i<map.size(); ++i)
        {
        	map[i].getDataProcessing().push_back(dp);          
      	}
      }
      
      ///Returns the the data processing information 
      DataProcessing getProcessingInfo_(DataProcessing::ProcessingAction action) const;

      ///Returns the the data processing information 
      DataProcessing getProcessingInfo_(const std::set<DataProcessing::ProcessingAction>& actions) const;

      //@}
      
			/// get DocumentIDTagger to assign DocumentIDs to maps
			const DocumentIDTagger& getDocumentIDTagger_() const;
			
			/**
				@brief Test mode 
			
				Test mode is enabled using the command line parameter @em -test .
				
				It disables writing of data, which would corrupt tests:
				- abolute paths (e.g. in consensus maps)
				- processing parameters (input/output files contain abolute paths as well)
				- current date
				- current OpenMS version
			*/
			bool test_mode_;

			/// .TOPP.ini file for storing system default parameters
			static String topp_ini_file_;
			
      /// width of console we are currently in (if not determinable, set to 80 as default)
      int console_width_;

      /// Debug level set by -debug
      Int debug_level_;

  };

} // namespace OpenMS

#endif //OPENMS_APPLICATIONS_TOPPBASE_H
