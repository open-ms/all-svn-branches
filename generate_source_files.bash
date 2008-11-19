#!/bin/bash

### This file generates source list files for CMake 
### each subdirectory which contains source file which
### need to be processed by CMake contains a seperate 
### file
### 11/2008 AB, CB

#### configs #########################################

SOURCE_FILENAME="sources.cmake"
SOURCE_FILE_EXTENSION="C"
HEADER_FILE_EXTENSION="h"
UIC_FILE_EXTENSION="ui"

######################################################


include_dirs=`ls -R include/OpenMS/* | grep ":" | cut -d ":" -f1`
source_dirs=`ls -R source/* | grep ":" | cut -d ":" -f1`
working_dir=`pwd`

### generate cmake files in include/
# MOC files
for d in ${include_dirs}
do
	source_name=${d}"/"${SOURCE_FILENAME}
	source_files=`cd ${d}; grep Q_OBJECT *.${HEADER_FILE_EXTENSION} 2> /dev/null | cut -d":" -f 1 2> /dev/null | sort | uniq; cd ${pwd}`
	if [ "${source_files}" == "" ] 
	then
		continue
	fi

  echo "Generating ${d}/${SOURCE_FILENAME}"

  # add directory name
  echo "### the directory name" > ${source_name}
  echo "set(directory ${d})" >> ${source_name}
  echo "" >> ${source_name}

  # add sources list
  echo "### list all filenames of the directory here" >> ${source_name}
  echo "set(sources_list" >> ${source_name}
  for f in ${source_files}
  do
    echo "${f}" >> ${source_name}
  done
  echo ")" >> ${source_name}
  echo "" >> ${source_name}

  echo "### add path to the filenames" >> ${source_name}
  echo "set(sources)" >> ${source_name}
  echo "foreach(i \${sources_list})" >> ${source_name}
  echo "  list(APPEND sources \${directory}/\${i})" >> ${source_name}
  echo "endforeach(i)" >> ${source_name}
  echo "" >> ${source_name}
	echo "### Apply MOC compiler" >> ${source_name}
	echo "QT4_WRAP_CPP(mocced_sources \${sources})" >> ${source_name}
	echo "" >> ${source_name}
  echo "### pass source file list to the upper instance" >> ${source_name}
  echo "set(OpenMS_sources \${OpenMS_sources} \${mocced_sources})" >> ${source_name}
  echo "" >> ${source_name}
done 

# UIC files
for d in ${include_dirs}
do
  source_name=${d}"/"${SOURCE_FILENAME}
  source_files=`cd ${d}; ls *.${UIC_FILE_EXTENSION} 2> /dev/null; cd ${pwd}`
  if [ "${source_files}" == "" ]
	then
    continue
  fi

  echo "Generating ${d}/${SOURCE_FILENAME}"

  # add directory name
  echo "### the directory name" > ${source_name}
  echo "set(directory ${d})" >> ${source_name}
  echo "" >> ${source_name}

  # add sources list
  echo "### list all filenames of the directory here" >> ${source_name}
  echo "set(sources_list" >> ${source_name}
  for f in ${source_files}
  do
    echo "${f}" >> ${source_name}
  done
  echo ")" >> ${source_name}
  echo "" >> ${source_name}

  echo "### add path to the filenames" >> ${source_name}
  echo "set(sources)" >> ${source_name}
  echo "foreach(i \${sources_list})" >> ${source_name}
  echo "  list(APPEND sources \${directory}/\${i})" >> ${source_name}
  echo "endforeach(i)" >> ${source_name}
  echo "" >> ${source_name}
  echo "### Apply UIC compiler" >> ${source_name}
  echo "QT4_WRAP_UI_OWN(uiced_sources \${sources})" >> ${source_name}
  echo "" >> ${source_name}
  echo "### pass source file list to the upper instance" >> ${source_name}
  echo "set(OpenMS_sources \${OpenMS_sources} \${uiced_sources})" >> ${source_name}
  echo "" >> ${source_name}
done

exit 1;

### generate cmake files in source/
for d in ${source_dirs}
do
	source_name=${d}"/"${SOURCE_FILENAME}
	source_files=`cd ${d}; ls *.${SOURCE_FILE_EXTENSION} 2> /dev/null; cd ${working_dir}`
	if [ "${source_files}" == "" ]
	then
		continue
	fi
	
	echo "Generating ${d}/${SOURCE_FILENAME}"

	# add directory name
	echo "### the directory name" > ${source_name}
	echo "set(directory ${d})" >> ${source_name}
	echo "" >> ${source_name}

	# add sources list
	echo "### list all filenames of the directory here" >> ${source_name}
	echo "set(sources_list" >> ${source_name}
	for f in ${source_files}
	do
		echo "${f}" >> ${source_name}
	done
	echo ")" >> ${source_name}
	echo "" >> ${source_name}

	echo "### add path to the filenames" >> ${source_name}
	echo "set(sources)" >> ${source_name}
	echo "foreach(i \${sources_list})" >> ${source_name}
  echo "	list(APPEND sources \${directory}/\${i})" >> ${source_name}
	echo "endforeach(i)" >> ${source_name}
	echo "" >> ${source_name}	
	echo "### pass source file list to the upper instance" >> ${source_name}
	echo "set(OpenMS_sources \${OpenMS_sources} \${sources})" >> ${source_name}
	echo "" >> ${source_name}
done

