// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/APPLICATIONS/TOPPToolParamHelper.h>

#include <QtGui/QGraphicsScene>
#include <QtGui/QMessageBox>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QRegExp>
#include <QtGui/QImage>

#include <QCoreApplication>

namespace OpenMS
{
	UInt TOPPASToolVertex::uid_ = 1;
	
	TOPPASToolVertex::TOPPASToolVertex()
		:	TOPPASVertex(),
      tool_name_(),
      tool_type_(),
			param_(),
			finished_(false),
			progress_color_(Qt::gray),
			iteration_nr_(0),
			input_list_length_(1)
		{
      pen_color_ = Qt::black;
      brush_color_ = QColor(245,245,245);
      if (TOPPToolParamHelper::initParam(param_, tool_name_, tool_type_, false, String())) {}
      connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
      connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
      connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
      connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
    }
	
	TOPPASToolVertex::TOPPASToolVertex(const String& name, const String& type, const String& tmp_path)
		: TOPPASVertex(),
      tool_name_(name),
      tool_type_(type),
			tmp_path_(tmp_path),
			param_(),
			finished_(false),
			progress_color_(Qt::gray),
			iteration_nr_(0),
			input_list_length_(1)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
    if (TOPPToolParamHelper::initParam(param_, tool_name_, tool_type_, false, String())) {}
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
	}
	
	TOPPASToolVertex::TOPPASToolVertex(const TOPPASToolVertex& rhs)
		:	TOPPASVertex(rhs),
      tool_name_(rhs.tool_name_),
      tool_type_(rhs.tool_type_),
			tmp_path_(rhs.tmp_path_),
			param_(rhs.param_),
			finished_(rhs.finished_),
			progress_color_(rhs.progress_color_),
			iteration_nr_(rhs.iteration_nr_),
			input_list_length_(rhs.input_list_length_)
	{
		pen_color_ = Qt::black;
		brush_color_ = QColor(245,245,245);
		connect (this, SIGNAL(toolStarted()), this, SLOT(toolStartedSlot()));
		connect (this, SIGNAL(toolFinished()), this, SLOT(toolFinishedSlot()));
		connect (this, SIGNAL(toolFailed()), this, SLOT(toolFailedSlot()));
		connect (this, SIGNAL(toolCrashed()), this, SLOT(toolCrashedSlot()));
	}

	TOPPASToolVertex::~TOPPASToolVertex()
	{
	
	}
	
	TOPPASToolVertex& TOPPASToolVertex::operator= (const TOPPASToolVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		
		param_ = rhs.param_;
    tool_name_ = rhs.tool_name_;
    tool_type_ = rhs.tool_type_;
		tmp_path_ = rhs.tmp_path_;
		finished_ = rhs.finished_;
		progress_color_ = rhs.progress_color_;
		iteration_nr_ = rhs.iteration_nr_;
		input_list_length_ = rhs.input_list_length_;
		
		return *this;
	}
	

	
	void TOPPASToolVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		editParam();
	}
	
	void TOPPASToolVertex::editParam()
	{
		QWidget* parent_widget = qobject_cast<QWidget*>(scene()->parent());
		String default_dir = "";
		
		// use a copy for editing
		Param edit_param(param_);
		
		QVector<Param::ParamEntry> hidden_entries;
		// remove type (should not be edited)
		if (edit_param.exists("type"))
		{
			hidden_entries.push_back(edit_param.getEntry("type"));
			edit_param.remove("type");
		}
		// remove entries that are handled by edges already, user should not see them
    QVector<TOPPIOInfo> input_infos;
		getInputParameters(input_infos);
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int index = (*it)->getTargetInParam();
			if (index < 0)
			{
				continue;
			}
			
			const String& name = input_infos[index].param_name;
			if (edit_param.exists(name))
			{
				hidden_entries.push_back(edit_param.getEntry(name));
				edit_param.remove(name);
			}
		}
		
    QVector<TOPPIOInfo> output_infos;
		getOutputParameters(output_infos);
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			int index = (*it)->getSourceOutParam();
			if (index < 0)
			{
				continue;
			}
			
			const String& name = output_infos[index].param_name;
			if (edit_param.exists(name))
			{
				hidden_entries.push_back(edit_param.getEntry(name));
				edit_param.remove(name);
			}
		}
		
    TOPPASToolConfigDialog dialog(parent_widget, edit_param, default_dir, tool_name_, tool_type_, hidden_entries);
		if (dialog.exec())
		{
			param_ = edit_param;
			// restore the removed entries
			foreach (const Param::ParamEntry& pe, hidden_entries)
			{
				StringList tags;
				for (std::set<String>::const_iterator it = pe.tags.begin(); it != pe.tags.end(); ++it)
				{
					tags.push_back(*it);
				}
				param_.setValue(pe.name, pe.value, pe.description, tags);
			}
			
			progress_color_ = Qt::gray;
			emit somethingHasChanged();
		}
		
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
	}
	
  void TOPPASToolVertex::getInputParameters(QVector<TOPPIOInfo>& input_infos)
	{
    TOPPToolParamHelper::getInputParameters(param_, input_infos);
	}
	
  void TOPPASToolVertex::getOutputParameters(QVector<TOPPIOInfo>& output_infos)
	{
    TOPPToolParamHelper::getOutputParameters(param_, output_infos);
	}
	
	void TOPPASToolVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		//painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform);
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
			pen.setColor(Qt::darkBlue);
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
		
		QPainterPath path;
		path.addRect(-70.0, -60.0, 140.0, 120.0);		
 		painter->drawPath(path);
 		
 		pen.setColor(pen_color_);
 		painter->setPen(pen);
    if (tool_type_ == "")
		{
      QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, tool_name_.toQString());
      painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), tool_name_.toQString());
		}
		else
		{
      QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, tool_name_.toQString());
      painter->drawText(-(int)(text_boundings.width()/2.0), -(int)(text_boundings.height()/3.0), tool_name_.toQString());
      text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, tool_type_.toQString());
      painter->drawText(-(int)(text_boundings.width()/2.0), +(int)(text_boundings.height()/1.33), tool_type_.toQString());
		}
		
		// progress light
		painter->setPen(Qt::black);
		painter->setBrush(progress_color_);
		painter->drawEllipse(46,-52, 14, 14);
		
		//topo sort number
		qreal x_pos = -64.0;
		qreal y_pos = -41.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
		
		if (progress_color_ != Qt::gray)
		{
			QString text;
			if (in_parameter_has_list_type_)
			{
				text = QString::number(iteration_nr_ == 1 ? input_list_length_ : 0);
				text += QString(" / ") + QString::number(input_list_length_); 
			}
			else
			{
				text = QString::number(iteration_nr_)+" / "+QString::number(num_iterations_);
			}
			QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
			painter->drawText((int)(62.0-text_boundings.width()), 48, text);
		}
	}
	
	QRectF TOPPASToolVertex::boundingRect() const
	{
		return QRectF(-71,-61,142,122);
	}
	
	QPainterPath TOPPASToolVertex::shape () const
	{
		QPainterPath shape;
		shape.addRect(-71.0, -61.0, 142.0, 122.0);				
		return shape;
	}
	
	const String& TOPPASToolVertex::getName()
	{
    return tool_name_;
	}
	
	const String& TOPPASToolVertex::getType()
	{
    return tool_type_;
	}
	
	void TOPPASToolVertex::runToolIfInputReady()
	{
		__DEBUG_BEGIN_METHOD__
		
		//check if everything ready
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (tv && !tv->isFinished())
			{
				// some tool that we depend on has not finished execution yet --> do not start yet
				debugOut_("Not run (parent not finished)");
				
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		// all inputs are ready --> GO!
		updateCurrentOutputFileNames();
		
		createDirs();
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QString ini_file = File::getTempDirectory().toQString()
							+QDir::separator()
							+getOutputDir().toQString()
							+QDir::separator()
              +tool_name_.toQString();
    if (tool_type_ != "")
		{
      ini_file += "_"+tool_type_.toQString();
		}
		ini_file += ".ini";
    TOPPToolParamHelper::writeParam(param_, tool_name_, ini_file);
		
		QStringList shared_args;
		shared_args	<< "-ini"
								<< ini_file
								<< "-no_progress";
    if (tool_type_ != "")
		{
      shared_args << "-type" << tool_type_.toQString();
		}
		
		ts->setPipelineRunning(true);
		emit toolStarted();
		
		/* Decide whether or not to iterate over the single files of the incoming lists
		 * depending on the type of the input parameter "-in" (file or list). */
		num_iterations_ = in_parameter_has_list_type_ ? 1 : input_list_length_;
		iteration_nr_ = 0; // needed in executionFinished()
		
		for (int i = 0; i < num_iterations_; ++i)
		{
			debugOut_(String("Enqueueing process nr ")+i);
			QStringList args = shared_args;
			
			// add all input file parameters
      QVector<TOPPIOInfo> in_params;
			getInputParameters(in_params);
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				int param_index = (*it)->getTargetInParam();
				if (param_index < 0)
				{
					std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
					break;
				}
				args << "-"+(in_params[param_index].param_name).toQString();
				
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
				if (tv)
				{
					int out_param_index = (*it)->getSourceOutParam();
					if (out_param_index < 0)
					{
						std::cerr << "TOPPAS: Output parameter index out of bounds!" << std::endl;
						break;
					}
					const QStringList& source_out_files = tv->current_output_files_[out_param_index];
					
					if (in_parameter_has_list_type_)
					{
						args << source_out_files;
					}
					else
					{
						if (i >= source_out_files.size())
						{
							std::cerr << "TOPPAS: Input list too short!" << std::endl;
							break;
						}
						args << source_out_files[i];
					}
					continue;
				}
				
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
				if (mv)
				{
					const QStringList& input_files = mv->getCurrentOutputList();

					if (in_parameter_has_list_type_)
					{
						args << input_files;
					}
					else
					{
						if (i >= input_files.size())
						{
							std::cerr << "TOPPAS: Input list too short!" << std::endl;
							break;
						}
						args << input_files[i];
					}
					continue;
				}

				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
					const QStringList& input_files = iflv->getFilenames();
					if (in_parameter_has_list_type_)
					{
						args << input_files;
					}
					else
					{
						if (i >= input_files.size())
						{
							std::cerr << "TOPPAS: Input list too short!" << std::endl;
							break;
						}
						args << input_files[i];
					}
					continue;
				}
			}
			
			// add all output file parameters
      QVector<TOPPIOInfo> out_params;
			getOutputParameters(out_params);
			
			for (int j = 0; j < out_params.size(); ++j)
			{
				// search for an out edge for this parameter
				for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
				{
					int param_index = (*it)->getSourceOutParam();
					
					if (j == param_index)
					{
						args << "-"+(out_params[param_index].param_name).toQString();
						const QStringList& output_files = current_output_files_[param_index];
						if (in_parameter_has_list_type_)
						{
							args << output_files;
						}
						else
						{
							if (i >= output_files.size())
							{
								std::cerr << "TOPPAS: Output list too short!" << std::endl;
								break;
							}
							args << output_files[i];
						}
						
						break; // (regardless of the number of out edges, every argument must appear only once)
					}
				}
			}
			
			//create process
			QProcess* p = new QProcess();
			p->setProcessChannelMode(QProcess::MergedChannels);
			connect(p,SIGNAL(finished(int,QProcess::ExitStatus)),this,SLOT(executionFinished(int,QProcess::ExitStatus)));
			connect(p,SIGNAL(readyReadStandardOutput()),this,SLOT(forwardTOPPOutput()));
			connect(ts,SIGNAL(terminateCurrentPipeline()),p,SLOT(kill()));
			
			//enqueue process
      std::cout << "TOPPAS: Enqueue: " << tool_name_ << " " << String(args.join(" ")) << std::endl;
      ts->enqueueProcess(p, tool_name_.toQString(), args);
		}
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASToolVertex::executionFinished(int ec, QProcess::ExitStatus es)
	{
		__DEBUG_BEGIN_METHOD__
		
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		if (es != QProcess::NormalExit)
		{
			ts->setPipelineRunning(false);
			emit toolCrashed();
			//clean up
			QProcess* p = qobject_cast<QProcess*>(QObject::sender());
			if (p)
			{
				delete p;
			}
			__DEBUG_END_METHOD__
			return;
		}
		
		if (ec != 0)
		{
			ts->setPipelineRunning(false);
			emit toolFailed();
			//clean up
			QProcess* p = qobject_cast<QProcess*>(QObject::sender());
			if (p)
			{
				delete p;
			}
			__DEBUG_END_METHOD__
			return;
		}
		
		++iteration_nr_;
		update(boundingRect());
		debugOut_(String("Increased iteration_nr_ to ")+iteration_nr_+" / "+num_iterations_);
		
		// notify the scene that this process has finished (so the next pending one can run)
		ts->runNextProcess();
		
		if (iteration_nr_ == num_iterations_) // all iterations performed --> proceed in pipeline
		{
			debugOut_("All iterations finished!");
			
			finished_ = true;
			emit toolFinished();
			
			if (all_written_output_files_.size() != current_output_files_.size())
			{
				all_written_output_files_.resize(current_output_files_.size());
			}
			for (int i = 0; i < current_output_files_.size(); ++i)
			{
				all_written_output_files_[i] << current_output_files_[i];
			}
			
			// notify all childs that we are finished, proceed in pipeline
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				TOPPASVertex* tv = (*it)->getTargetVertex();
				debugOut_(String("Starting child ")+tv->getTopoNr());
				
				TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
				if (ttv)
				{
					ttv->runToolIfInputReady();
					continue;
				}
				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(tv);
				if (mv)
				{
					mv->forwardPipelineExecution();
					continue;
				}
				TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
				if (oflv)
				{
					oflv->finish();
					continue;
				}
			}
			
			debugOut_("All children started!");
		}
		
		//clean up
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (p)
		{
			delete p;
		}
		
		__DEBUG_END_METHOD__
	}
	
	bool TOPPASToolVertex::isFinished()
	{
		return finished_;
	}
	
	const Param& TOPPASToolVertex::getParam()
	{
		return param_;
	}
	
	void TOPPASToolVertex::setParam(const Param& param)
	{
		param_ = param;
	}
	
	const QVector<QStringList>& TOPPASToolVertex::getCurrentOutputFileNames()
	{
		return current_output_files_;
	}
	
	const QVector<QStringList>& TOPPASToolVertex::getAllWrittenOutputFileNames()
	{
		return all_written_output_files_;
	}
	
	void TOPPASToolVertex::updateCurrentOutputFileNames()
	{
    QVector<TOPPIOInfo> in_params;
		input_list_length_ = 1; // stays like that if -in param is not a list
		bool found_in_parameter = false;
		getInputParameters(in_params);
		QStringList input_file_basenames;
		
		bool force = false;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			int param_index = (*it)->getTargetInParam();
			if (param_index < 0)
			{
				std::cerr << "TOPPAS: Input parameter index out of bounds!" << std::endl;
				break;
			}
			
			if (in_params[param_index].param_name == "in" || force)
			{
				found_in_parameter = true;
        in_parameter_has_list_type_ = (in_params[param_index].type == TOPPIOInfo::IOT_LIST);
				
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>((*it)->getSourceVertex());
				if (iflv)
				{
					const QStringList& input_files = iflv->getFilenames();
					input_list_length_ = input_files.count();
					foreach (const QString& str, input_files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
					break;
				}

				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
				if (tv)
				{
					int out_param_index = (*it)->getSourceOutParam();
					if (out_param_index < 0)
					{
						std::cerr << "TOPPAS: Output parameter index out of bounds!" << std::endl;
						break;
					}
					const QStringList& input_files = tv->current_output_files_[out_param_index];
					input_list_length_ = input_files.count();
					foreach (const QString& str, input_files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
					break;
				}

				TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
				if (mv)
				{
					const QStringList& files = mv->getCurrentOutputList();
					input_list_length_ = files.count();

					foreach (const QString& str, files)
					{
						input_file_basenames.push_back(File::basename(str).toQString());
					}
					break;
				}
			}
			
			//if last iteration and still no "in" parameter found, repeat last iteration and treat the edge as input parameter (dirty - TODO)
			if (it == inEdgesEnd()-1 && !found_in_parameter)
			{
				--it;
				force = true;
			}
		}
		
		// now, update the output file names:
    QVector<TOPPIOInfo> out_params;
		getOutputParameters(out_params);
		
		current_output_files_.clear();
		current_output_files_.resize(out_params.size());
		
		for (int i = 0; i < out_params.size(); ++i)
		{
			// search for an out edge for this parameter
			for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
			{
				int param_index = (*it)->getSourceOutParam();
				
				if (i == param_index) // corresponding out edge found
				{
					// check if tool consumes list and outputs single file (such as IDMerger or FileMerger)
          if (in_parameter_has_list_type_ && out_params[param_index].type == TOPPIOInfo::IOT_FILE)
					{
						QString f = File::getTempDirectory().toQString()
							+QDir::separator()
							+getOutputDir().toQString()
							+QDir::separator()
							+out_params[param_index].param_name.toQString()
							+QDir::separator()
							+input_file_basenames.first()
							+"_to_"
							+input_file_basenames.last()
							+"_merged_tmp"+QString::number(uid_++);
						current_output_files_[param_index].push_back(f);
					}
					else
					{
						foreach (const QString& str, input_file_basenames)
						{
							QString f = File::getTempDirectory().toQString()
								+QDir::separator()
								+getOutputDir().toQString()
								+QDir::separator()
								+out_params[param_index].param_name.toQString()
								+QDir::separator()
								+str;
							QRegExp rx("_tmp\\d+$");
							int tmp_index = rx.indexIn(f);
							//std::cout << "tmp_index: " << tmp_index << std::endl;
							if (tmp_index != -1)
							{
								f = f.left(tmp_index);
							}
							f += "_tmp" + QString::number(uid_++);
							current_output_files_[param_index].push_back(f);
						}
					}

					break; // we are done, don't push_back further file names for this parameter
				}
			}
		}
	}

	void TOPPASToolVertex::forwardTOPPOutput()
	{
		QProcess* p = qobject_cast<QProcess*>(QObject::sender());
		if (!p)
		{
			return;
		}
		
		QString out = p->readAllStandardOutput();
		emit toppOutputReady(out);
	}
	
	void TOPPASToolVertex::setProgressColor(const QColor& c)
	{
		progress_color_ = c;
	}
	
	QColor TOPPASToolVertex::getProgressColor()
	{
		return progress_color_;
	}
	
	void TOPPASToolVertex::toolStartedSlot()
	{
		progress_color_ = Qt::yellow;
		update(boundingRect());
	}
	
	void TOPPASToolVertex::toolFinishedSlot()
	{
		progress_color_ = Qt::green;
		update(boundingRect());
	}
	
	void TOPPASToolVertex::toolFailedSlot()
	{
		progress_color_ = Qt::red;
		update(boundingRect());
	}

	void TOPPASToolVertex::toolCrashedSlot()
	{
		progress_color_ = Qt::red;
		update(boundingRect());
	}

	void TOPPASToolVertex::inEdgeHasChanged()
	{
    // something has changed --> tmp files might be in --> reset
		reset(true);
		
		TOPPASVertex::inEdgeHasChanged();
	}
	
	void TOPPASToolVertex::openInTOPPView()
	{
    QVector<TOPPIOInfo> out_infos;
		getOutputParameters(out_infos);
		
		if (out_infos.size() == all_written_output_files_.size())
		{
			foreach (const QStringList& files, all_written_output_files_)
			{
				if (files.size() > 0)
				{
					QProcess* p = new QProcess();
					p->setProcessChannelMode(QProcess::ForwardedChannels);
          QString toppview_executable;
          toppview_executable = "TOPPView";

          p->start(toppview_executable, files);
          if(!p->waitForStarted())
          {
            // execution failed
            std::cerr << p->errorString().toStdString() << std::endl;
#if defined(Q_WS_MAC)
            std::cerr << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
#endif

          }
				}
			}
		}
	}
	
	String TOPPASToolVertex::getOutputDir()
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		String workflow_dir = File::removeExtension(File::basename(ts->getSaveFileName()));
		if (workflow_dir == "")
		{
			workflow_dir = "Untitled_workflow";
		}
		String dir = String("TOPPAS_tmp")+
			String(QDir::separator())+
			workflow_dir+
			String(QDir::separator())+
			get3CharsNumber_(topo_nr_)+"_"+getName();
		if (getType() != "")
		{
			dir += "_"+getType();
		}
		
		return dir;
	}
	
	void TOPPASToolVertex::createDirs()
	{
		QDir current_dir(File::getTempDirectory().toQString());
		
		if (!current_dir.mkpath(getOutputDir().toQString()))
		{
			std::cerr << "TOPPAS: Could not create path " << getOutputDir() << std::endl;
		}
		
		foreach (const QStringList& files, current_output_files_)
		{
			if (!files.isEmpty())
			{
				QString dir = File::path(files.first()).toQString();
				if (!File::exists(dir))
				{
					if (!current_dir.mkpath(dir))
					{
						std::cerr << "TOPPAS: Could not create path " << String(dir) << std::endl;
					}
				}
			}
		}
	}
	
	void TOPPASToolVertex::setTopoNr(UInt nr)
	{
		if (topo_nr_ != nr)
		{
			// topological number changes --> output dir changes --> reset
			reset(true);
			topo_nr_ = nr;
			emit somethingHasChanged();
		}
	}
	
	void TOPPASToolVertex::reset(bool reset_all_files)
	{
		__DEBUG_BEGIN_METHOD__
		
		finished_ = false;
		current_output_files_.clear();
		progress_color_ = Qt::gray;
		
		if (reset_all_files)
		{
			all_written_output_files_.clear();
			QString remove_dir = File::getTempDirectory().toQString() + QDir::separator() + getOutputDir().toQString();
			if (File::exists(remove_dir))
			{
				removeDirRecursively_(remove_dir);
			}
			// reset UID for tmp files
			uid_ = 1;
		}
		
		TOPPASVertex::reset(reset_all_files);
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASToolVertex::checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run)
	{
		__DEBUG_BEGIN_METHOD__
		
		//all parents checked?
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			if (!((*it)->getSourceVertex()->isScListLengthChecked()))
			{
				__DEBUG_END_METHOD__
				return;
			}
		}
		
		// do all input lists have equal length?
		int parent_per_round = (*inEdgesBegin())->getSourceVertex()->getScFilesPerRound();
		int parent_total = (*inEdgesBegin())->getSourceVertex()->getScFilesTotal();		
		
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			if ((*it)->getSourceVertex()->getScFilesPerRound() != parent_per_round)
			{
				unequal_per_round.push_back(QString::number(topo_nr_));
				break;
			}
		}
		
		// n:1 tool? --> files per round = 1
    QVector<TOPPIOInfo> input_infos;
		getInputParameters(input_infos);
    QVector<TOPPIOInfo> output_infos;
		getOutputParameters(output_infos);
		bool in_param_list_type = false;
		bool out_param_file_type = false;
    foreach (const TOPPIOInfo& io, input_infos)
		{
      if (io.param_name == "in" && io.type == TOPPIOInfo::IOT_LIST)
			{
				in_param_list_type = true;
			}
		}
    foreach (const TOPPIOInfo& io, output_infos)
		{
      if (io.param_name == "out" && io.type == TOPPIOInfo::IOT_FILE)
			{
				out_param_file_type = true;
			}
		}
		
		if (in_param_list_type && out_param_file_type)
		{
			sc_files_per_round_ = 1;
			sc_files_total_ = parent_total / parent_per_round;
		}
		else
		{
			sc_files_per_round_ = parent_per_round;
			sc_files_total_ = parent_total;
		}
		
		sc_list_length_checked_ = true;
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			tv->checkListLengths(unequal_per_round, unequal_over_entire_run);
		}
		
		__DEBUG_END_METHOD__
	}

  bool TOPPASToolVertex::refreshParameters()
  {
    return TOPPToolParamHelper::refreshParameters(param_, tool_name_, tool_type_);
  }
}

