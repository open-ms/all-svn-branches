#include <OpenMS/APPLICATIONS/TOPPToolParamHelper.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QMessageBox>
#include <QtCore/QString>
#include <QtCore/QDir>
#include <QtCore/QFile>

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

    bool changed = initParam(tool_parameter, tool_name, tool_type, old_ini_file);

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

  bool TOPPToolParamHelper::initParam(Param& tool_param, String tool_name, String tool_type, const String& old_ini_file)
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
        QMessageBox::critical(0,"Error",(String("Could not open '")+old_ini_file+"'!").c_str());
        return false;
      }
      call += " -ini " + String(old_ini_file);
    }

    if (system(call.c_str()) != 0)
    {
      QMessageBox::critical(0,"Error",(String("Could not execute '")+call+"'!\n\nMake sure the TOPP tools are in your $PATH variable, that you have write permission in the temporary file path, and that there is space left in the temporary file path.").c_str());
      return false;
    }
    if(!File::exists(ini_file))
    {
      QMessageBox::critical(0,"Error",(String("Could not open '")+ini_file+"'!").c_str());
      return false;
    }

    tmp_param.load(String(ini_file).c_str());
    tool_param = tmp_param.copy(tool_name + ":1:",true);
    //param_.remove("log");
    //param_.remove("no_progress");
    //param_.remove("debug");
    //// handled by TOPPAS anyway:
    //param_.remove("type");

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
        StringList valid_types;

        const String& desc = it->description;
        String::SizeType index = desc.find("valid formats",0);
        if (index != String::npos)
        {
          String::SizeType types_start_pos = desc.find("'",index) + 1;
          String::SizeType types_length = desc.find("'",types_start_pos) - types_start_pos;
          String types_string = desc.substr(types_start_pos, types_length);
          if (types_string.find(",",0) == String::npos)
          {
            valid_types.push_back(types_string.trim());
          }
          else
          {
            types_string.split(',', valid_types);
          }
        }

        TOPPIOInfo io_info;
        io_info.param_name = it->name;
        io_info.valid_types = valid_types;
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
}
