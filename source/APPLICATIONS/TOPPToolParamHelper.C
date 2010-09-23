#include <OpenMS/APPLICATIONS/TOPPToolParamHelper.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtGui/QMessageBox>
#include <QtCore/QString>
#include <QtCore/QDir>
#include <QtCore/QFile>

using namespace std;

namespace OpenMS
{

  bool TOPPToolParamHelper::refreshParameters(Param& tool_parameter, const String& tool_name, const String& tool_type)
  {
    QString old_ini_file = QDir::tempPath() + QDir::separator() + "TOPPAS_" + tool_name.toQString() + "_";
    if (tool_type != "")
    {
      old_ini_file += tool_type.toQString() + "_";
    }
    old_ini_file += File::getUniqueName().toQString() + "_tmp_OLD.ini";
    writeParam(tool_parameter, tool_name, old_ini_file);

    bool changed = initParam(tool_parameter, tool_name, tool_type, false, old_ini_file);

    return changed;
  }

  void TOPPToolParamHelper::writeParam(const Param& param, const String& tool_name, const QString& ini_file)
  {
    Param save_param;
    save_param.setValue(tool_name+":1:toppas_dummy", DataValue("blub"));
    save_param.insert(tool_name+":1:", param);
    save_param.remove(tool_name+":1:toppas_dummy");
    save_param.setSectionDescription(tool_name+":1", "Instance '1' section for '"+tool_name+"'");
    save_param.store(ini_file);
  }

  StringList TOPPToolParamHelper::getToolTypes(String tool_name)
  {
    Map<String,StringList> tool_name_2_types = TOPPBase::getToolList();
    Map<String,StringList>::iterator it = tool_name_2_types.find(tool_name);
    if (it == tool_name_2_types.end())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Requesting TOPP tool types of unknown TOPP tool "+tool_name);
    }
    return it->second;
  }

  bool TOPPToolParamHelper::initParam(Param& tool_param, String tool_name, String tool_type, bool show_messagebox_on_error, const String& old_ini_file)
  {
    Param tmp_param;
    QString ini_file = QDir::tempPath() + QDir::separator() + "TOPPAS_" + tool_name.toQString() + "_";
    if (tool_type != "")
    {
      ini_file += tool_type.toQString() + "_";
    }
    ini_file += File::getUniqueName().toQString() + "_tmp.ini";

    String call = tool_name + " -write_ini " + ini_file;
    if (tool_type != "")
    {
      call += " -type " + tool_type;
    }
    if (old_ini_file != "")
    {
      if (!File::exists(old_ini_file))
      {
        if (show_messagebox_on_error)
        {
          QMessageBox::critical(0,"Error",(String("Could not open '")+old_ini_file+"'!").c_str());
        }
        return false;
      }
      call += " -ini " + String(old_ini_file);
    }

    if (system(call.c_str()) != 0)
    {
      if (show_messagebox_on_error)
      {
        QMessageBox::critical(0,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
      }
      return false;
    }
    if(!File::exists(ini_file))
    {
      if (show_messagebox_on_error)
      {
        QMessageBox::critical(0,"Error",(String("Could not open '")+ini_file+"'!").c_str());
      }
      return false;
    }

    tmp_param.load(String(ini_file).c_str());
    tool_param = tmp_param.copy(tool_name + ":1:",true);

    writeParam(tool_param, tool_name, ini_file);
    bool changed = false;
    if (old_ini_file != "")
    {
      //check if ini file has changed (quick & dirty by file size)
      QFile q_ini(ini_file);
      QFile q_old_ini(old_ini_file.toQString());
      changed = q_ini.size() != q_old_ini.size();
      QFile::remove(old_ini_file.toQString());
    }
    QFile::remove(ini_file);

    return changed;
  }

  void TOPPToolParamHelper::getInputParameters(const Param& param, QVector<TOPPIOInfo>& io_infos)
  {
    getParameters_(param, io_infos, true);
  }

  void TOPPToolParamHelper::getOutputParameters(const Param& param, QVector<TOPPIOInfo>& io_infos)
  {
    getParameters_(param, io_infos, false);
  }

  void TOPPToolParamHelper::getParameters_(const Param& param_, QVector<TOPPIOInfo>& io_infos, bool input_params)
  {
    String search_tag = input_params ? "input file" : "output file";

    io_infos.clear();

    for (Param::ParamIterator it = param_.begin(); it != param_.end(); ++it)
    {

      if (it->tags.count(search_tag))
      {
        StringList valid_file_types;

        const String& desc = it->description;
        String::SizeType index = desc.find("valid formats",0);
        if (index != String::npos)
        {
          String::SizeType types_start_pos = desc.find("'",index) + 1;
          String::SizeType types_length = desc.find("'",types_start_pos) - types_start_pos;
          String types_string = desc.substr(types_start_pos, types_length);
          if (types_string.find(",",0) == String::npos)
          {
            valid_file_types.push_back(types_string.trim());
          }
          else
          {
            types_string.split(',', valid_file_types);
          }
        }

        TOPPIOInfo io_info;
        io_info.param_name = it->name;
        io_info.valid_types = valid_file_types;
        if (it->value.valueType() == DataValue::STRING_LIST)
        {
          io_info.type = TOPPIOInfo::IOT_LIST;
        }
        else if (it->value.valueType() == DataValue::STRING_VALUE)
        {
          io_info.type = TOPPIOInfo::IOT_FILE;
        }
        else
        {
          std::cerr << "TOPPAS: Unexpected parameter value!" << std::endl;
        }
        io_infos.push_back(io_info);
      }
    }
    // order in param can change --> sort
    qSort(io_infos);
  }

  bool TOPPToolParamHelper::toolAcceptsFileExtension(const Param& target_tool_params, String extension)
  {
    // extract input parameters of target TOPP tool from param
    QVector<TOPPIOInfo> io_infos;
    TOPPToolParamHelper::getInputParameters(target_tool_params, io_infos);

    #ifdef DEBUG_TOPPVIEW
    cout << "Number of Input slots: " << io_infos.size() << endl;
    #endif

    // print input slots
    Int in_index = -1;
    for(int i=0; i!=io_infos.size(); ++i)
    {

      #ifdef DEBUG_TOPPVIEW
      cout << "in slot " << i << " " << io_infos[i].param_name << endl;
      cout << "in types: " << io_infos[i].valid_types << endl;
      #endif

      if (io_infos[i].param_name == "in" && io_infos[i].type == TOPPIOInfo::IOT_FILE)
      {
        in_index = i;
      }
    }

    #ifdef DEBUG_TOPPVIEW
    cout << "-in is on slot: " << in_index << endl;
    #endif

    // no proper in parameter found that accepts a single file
    if (in_index == -1)
    {
      return false;
    }

    StringList target_param_types = io_infos[in_index].valid_types;

    if (target_param_types.empty())
    {
      return false;
    }

    extension.toLower();

    #ifdef DEBUG_TOPPVIEW
    cout << "ext: " << extension << endl;
    #endif

    // check file type compatibility using the extension
    bool mismatch = true;
    for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
    {
      String other_ext = *it;
      other_ext.toLower();
      #ifdef DEBUG_TOPPVIEW
      cout << "other_ext: " << other_ext << endl;
      #endif

      if (extension == other_ext || extension == "gz")
      {
        mismatch = false;
        break;
      }
    }

    if (mismatch)
    {
      return false;
    }
    return true;
  }


  bool TOPPToolParamHelper::toolDeliversFileExtension(const Param& target_tool_params, StringList supported_outfile_extensions)
  {
    // extract input parameters of target TOPP tool from param
    QVector<TOPPIOInfo> io_infos;
    TOPPToolParamHelper::getOutputParameters(target_tool_params, io_infos);
    cout << "ios size out: " << io_infos.size() << endl;

    // tool has no output?
    if (io_infos.size()==0)
    {
      return false;
    }

    // find in parameter of possibly several in parameters
    Int out_index = -1;
    for(int i=0; i!=io_infos.size(); ++i)
    {
      #ifdef DEBUG_TOPPVIEW
      cout << "out slot " << i << " " << io_infos[i].param_name << endl;
      cout << "valid types: " << io_infos[i].valid_types << endl;
      #endif

      if (io_infos[i].param_name == "out")
      {
        out_index = i;
      }
    }

    // out slot not found?
    if (out_index == -1)
    {
      return false;
    }

    // get valid types of "out" parameter
    StringList target_param_types = io_infos[out_index].valid_types;

    if (target_param_types.empty())
    {
      return false;
    }

    for(UInt i=0; i!= supported_outfile_extensions.size(); ++i)
    {
      supported_outfile_extensions[i].toLower();
    }

    // check file type compatibility using the extension
    bool mismatch = true;
    for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
    {
      String other_ext = *it;
      other_ext.toLower();

      #ifdef DEBUG_TOPPVIEW
      cout << "out other_ext: " << other_ext << endl;
      #endif

      if (std::find(supported_outfile_extensions.begin(), supported_outfile_extensions.end(), other_ext) != supported_outfile_extensions.end())
      {
        mismatch = false;
        break;
      }
    }

    if (mismatch)
    {
      return false;
    }
    return true;
  }

}
