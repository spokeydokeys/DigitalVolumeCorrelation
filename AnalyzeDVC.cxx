//      AnalyzeDVC.cxx
//      
//      Copyright 2012 Seth Gilchrist <seth@seth-desktop-CHHM6>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.
//      
//      

#ifndef ANALYZEDVC_H
#define ANALYZEDVC_H

#include <iostream>
#include "DICMesh.cxx"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

// the following provides updates for the RegularStep optimizer
class CommandIterationUpdate : public itk::Command
{
public:
typedef  CommandIterationUpdate   Self;
typedef  itk::Command             Superclass;
typedef itk::SmartPointer<Self>  Pointer;
itkNewMacro( Self );

protected:
CommandIterationUpdate() {};

public:

typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
typedef const OptimizerType                         *OptimizerPointer;
std::string											m_LogfileName;

void Execute(itk::Object *caller, const itk::EventObject & event)
  {
	Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
	OptimizerPointer optimizer = 
						 dynamic_cast< OptimizerPointer >( object );

	if( ! itk::IterationEvent().CheckEvent( &event ) )
	  {
	  return;
	  }
	
	std::stringstream msg("");
	msg << optimizer->GetCurrentIteration() << " = "  << optimizer->GetValue() << " : " << optimizer->GetCurrentPosition()<<" step size "<<optimizer->GetCurrentStepLength() << std::endl;
	this->WriteToLogfile( msg.str() );
}

void SetLogfileName( std::string logfileName )
{
	this->m_LogfileName = logfileName;
}

void WriteToLogfile( std::string characters )
{
	std::ofstream outFile;
	outFile.open(this->m_LogfileName.c_str(), std::ofstream::app);
	if(!outFile.is_open())
	{
		std::cerr<<"Logfile error!  Cannot open file."<<std::endl;
		std::abort();
	}
	std::cout<< characters;
	outFile << characters;
	
	outFile.close();
}   
};


template<typename TFixedImage, typename TMovingImage>
class AnalyzeDVC : public DICMesh<TFixedImage, TMovingImage>
{
public:

typedef double					StepLengthType;
typedef std::string				ConfigurationFileNameType;
typedef TFixedImage				FixedImageType;
typedef TMovingImage			MovingImageType;

/** Constructor  **/
AnalyzeDVC()
{
	this->GetRegistrationMethod()->SetNumberOfThreads( 2 );	// default to 2 threads
	m_configFileName.clear(); 							// must be set by user
	
	m_GlobalMaxStep = 0.010;						// must be set by user
	m_GlobalMinStep = 0.005;						// must be set by user
	
	m_InitialDVCMaxStep = 0.010;					// must be set by user
	m_InitialDVCMinStep = 0.0005;					// must be set by user
	m_SecondaryDVCMaxStep = 0;						// must be set by user
	m_SecondaryDVCMinStep = 0;						// must be set by user
	//~ m_TertiaryDVCMaxStep = 0;						// must be set by user
	//~ m_TertiaryDVCMinStep = 0;						// must be set by user
	
	m_IdispErrorToll = 2;							// default to 2 stdev
	m_IdispReplaceSigma = 0.5;						// must be set by user
	m_IdispReplaceMean = 0;							// must be set by user
	m_IdispSmoothSigma = 0.1;						// must be set by user
	m_IdispSmoothMean = 0;							// must be set by user
	m_IstrainErrorToll = 1;							// default to 1 stdev
	m_IstrainReplaceSigma = 0;						// must be set by user
	m_IstrainReplaceMean = 0;						// must be set by user
	m_IstrainSmoothSigma = 0;						// must be set by user
	m_IstrainSmoothMean = 0;							// must be set by user

	m_SdispErrorToll = 2;							// default to 2 stdev
	m_SdispReplaceSigma = 0;						// must be set by user
	m_SdispReplaceMean = 0;							// must be set by user
	m_SdispSmoothSigma = 0;							// must be set by user
	m_SdispSmoothMean = 0;							// must be set by user
	m_SstrainErrorToll = 1;							// default to 1 stdev
	m_SstrainReplaceSigma = 4;						// must be set by user
	m_SstrainReplaceMean = 0;						// must be set by user
	m_SstrainSmoothSigma = 2;						// must be set by user
	m_SstrainSmoothMean = 0;						// must be set by user
		
	m_SecondaryDVC = 0;								// default to forgo secondary DVC
	//~ m_TertiaryDVC = 0;								// default to forgo tertiary DVC
	
	m_RestartFile = 0;								// default to not use the restart methods
	
	m_fixedFileName.clear();						// must be set by user
	m_movingFileName.clear();						// must be set by user
	m_meshFileName.clear();							// must be set by user
	m_outputDirectory.clear();						// must be set by user
	
	m_observer = CommandIterationUpdate::New();
}

/** Deconstructor **/
~AnalyzeDVC() {}

/** Set Configuration File **/
void SetConfigurationFile(std::string fileName )
{
	this->m_configFileName = fileName;
}

/** A function to read the configuration file.  Returns true if
 * the file is read successfully. False if there is an error. **/
 /* Config File Format is a key = value style. Valid keys and value types are below. 
  * Line endings signifiy the end of a value.
  * Lines starting with '#' are comments."
  * Values shown below in () are defaults.

# Fixed image file name
FIXEDIMAGEFILE=string (0)
# Moving image file name
MOVINGIMAGEFILE=string (0)
# Mesh image (gmsh or vtk) file name
MESHFILENAME=string (0)
# Output folder
OUTPUTFOLDER=string (0)
# Interrogation region radius
IRRADIUS=int (0)
# Number of threads
NTHREADS=int (2) 
# Max/Min step length for the global registration
GLOBALMAXSTEP=double (0.010)
GLOBALMINSTEP=double (0.005)
# Max/Min step length for initial DVC
INITIALDVCMAXSTEP=double (0.010)
INITIALDVCMINSTEP=double (0.0005)
# Max/Min step length for second DVC (if executing)
SECONDARYDVCMAXSTEP=double (0)
SECONDARYDVCMINSTEP=double (0)
# Flag to perform second DVC
PERFORMSECONDARYDVC=bool (0)
# Error detection and handeling after initial DVC
# Displacement error tollerance in stdev from neighbourhood mean
IDISPLACEMENTERRORTOLLERANCE=double (2)
# Displacement Replacement
IDISPREPLACESIGMA=double (0.5)
IDISPREPLACEMEAN=double (0)
# Displcement Smoothing
IDISPLACESMOOTHSIGMA=double (0.1)
IDISPLACESMOOTHMEAN=double (0)
# Strain error tolleranc in stdev from neighbourhood mean
ISTRAINERRORTOLLERANCE=double (1)
# Strain Replacement
ISTRAINREPLACESIGMA=double (0)
ISTRAINREPLACEMEAN=double (0)
#Strain smoothing
ISTRAINSMOOTHSIGMA=double (0)
ISTRAINSMOOTHMEAN=double (0)
# Error detection and handelling after secondary DVC
# Displacement error tollerance in stdev from neighbourhood mean
SDISPLACEMENTERRORTOLLERANCE=double (2)
# Displacement Replacement
SDISPREPLACESIGMA=double (0)
SDISPREPLACEMEAN=double (0)
# Displcement Smoothing
SDISPLACESMOOTHSIGMA=double (0)
SDISPLACESMOOTHMEAN=double (0)
# Strain error tolleranc in stdev from neighbourhood mean
SSTRAINERRORTOLLERANCE=double (1)
# Strain Replacement
SSTRAINREPLACESIGMA=double (4)
SSTRAINREPLACEMEAN=double (0)
# Strain smoothing
SSTRAINSMOOTHSIGMA=double (2)
SSTRAINSMOOTHMEAN=double (0)
* 
* (Currently not implementd)
* # Max/Min stop length for tertiary DVC 
* TERTIARYDVCMAXSTEP=double (0)
* TERTIARYDVCMAXSTEP=double (0)
* # Flag to perform tertiary DVC
* PERFORMTERTIARTYDVC=bool (0)
  * */
bool ReadConfigureationFile()
{
	if ( this->m_configFileName.empty() ){
		std::cerr<<	"Configuration file not specified."<<std::endl;
		return 1;
	}
	
	std::ifstream configFileInput( this->m_configFileName.c_str() ); // open file for reading
	if(!configFileInput){ // if the file fails to open, give an error and abort
		std::cerr<<	"Cannot open configuration file for reading."<<std::endl <<
					"Please check the file name and permissions and try again."<<std::endl;
		return 1;
	}

	while ( !configFileInput.eof() ){
		std::string cLine;
		std::string key;
		std::string value;
		std::getline(configFileInput,cLine);
		
		// skip comment lines
		key = "#";
		if ( !cLine.compare(0,key.size(),key) ) {continue;}
		
		// if the fixed image name  TEST THIS OUT!!!!!
		key = "FIXEDIMAGEFILE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_fixedFileName = value;
			continue;
		}
		// if the moving image name
		key = "MOVINGIMAGEFILE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_movingFileName = value;
			continue;
		}
		// if mesh file name
		key = "MESHFILENAME";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_meshFileName = value;
			continue;
		}
		// if output folder
		key = "OUTPUTFOLDER";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_outputDirectory = value;
			this->SetOuputDirectory( this->m_outputDirectory );
			std::string logfile = this->m_outputDirectory+"/logfile.txt";
			this->SetLogfileName( logfile );
			continue;
		}
		
		// if IR radius
		key = "IRRADIUS";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->SetInterrogationRegionRadius( atoi( value.c_str()) );
			continue;
		}
		// if number of threads
		key = "NTHREADS";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->GetRegistrationMethod()->SetNumberOfThreads( atoi( value.c_str()) );
			continue;
		}
		
		// if Global max step size
		key = "GLOBALMAXSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_GlobalMaxStep = atof( value.c_str() );
			continue;
		}
		// if Global min step size
		key = "GLOBALMINSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_GlobalMinStep = atof( value.c_str() );
			continue;
		}
		
		// if initial max step size
		key = "INITIALDVCMAXSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_InitialDVCMaxStep = atof( value.c_str() );
			continue;
		}
		// if initial min step size
		key = "INITIALDVCMINSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_InitialDVCMinStep = atof( value.c_str() );
			continue;
		}
		
		// if secondary max step size
		key = "SECONDARYDVCMAXSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SecondaryDVCMaxStep = atof( value.c_str() );
			continue;
		}
		// if secondary min step size
		key = "SECONDARYDVCMINSTEP";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SecondaryDVCMinStep = atof(value.c_str());
			continue;
		}
		
		//~ // if tertiary max step size
		//~ key = "TERTIARYDVCMAXSTEP";
		//~ if ( !cLine.compare(0,key.size(),key) ){
			//~ value.assign(cLine,key.size()+1,511);
			//~ this->m_TertiaryDVCMaxStep = atof( value.c_str() );
			//~ continue;
		//~ }
		//~ // if tertiary min step size
		//~ key = "TERTIARYDVCMAXSTEP";
		//~ if ( !cLine.compare(0,key.size(),key) ){
			//~ value.assign(cLine,key.size()+1,511);
			//~ this->m_TertiaryDVCMinStep = atof( value.c_str() );
			//~ continue;
		//~ }
		
		// if secondary flag
		key = "PERFORMSECONDARYDVC";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SecondaryDVC = atoi( value.c_str() );
			continue;
		}
		//~ // if tertiary flag
		//~ key = "PERFORMTERTIARTYDVC";
		//~ if ( !cLine.compare(0,key.size(),key) ){
			//~ value.assign(cLine,key.size()+1,511);
			//~ this->m_TertiaryDVC = atoi( value.c_str() );
			//~ continue;
		//~ }
		
		// if intial displacement error tollerance
		key = "IDISPLACEMENTERRORTOLLERANCE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IdispErrorToll = atof( value.c_str() );
			continue;
		}
		// if after initial displacement replace sigma
		key = "IDISPREPLACESIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IdispReplaceSigma = atof( value.c_str() );
			continue;
		}
		// if after initial displacemnet replace mean
		key = "IDISPREPLACEMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IdispReplaceMean = atof( value.c_str() );
			continue;
		}
		
		// if after initial displacement smooth sigma
		key = "IDISPLACESMOOTHSIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IdispSmoothSigma = atof( value.c_str() );
			continue;
		}
		// if after intial displacement smooth mean
		key = "IDISPLACESMOOTHMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IdispSmoothMean = atof( value.c_str() );
			continue;
		}
		
		// if strain intial error tollerance
		key = "ISTRAINERRORTOLLERANCE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IstrainErrorToll = atof( value.c_str() );
			continue;
		}				
		// if after initial strain replace sigma
		key = "ISTRAINREPLACESIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IstrainReplaceSigma = atof( value.c_str() );
			continue;
		}
		// if after initial strain replace mean
		key = "ISTRAINREPLACEMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IstrainReplaceMean = atof( value.c_str() );
			continue;
		}
		// if after initil strain smooth sigma
		key = "ISTRAINSMOOTHSIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IstrainSmoothSigma = atof( value.c_str() );
			continue;
		}
		// if after initia strain replace mean
		key = "ISTRAINSMOOTHMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_IstrainSmoothMean = atof( value.c_str() );
			continue;
		}
		
		
		// if secondary displacement error tollerance
		key = "SDISPLACEMENTERRORTOLLERANCE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SdispErrorToll = atof( value.c_str() );
			continue;
		}
		// if after secondary displacement replace sigma
		key = "SDISPREPLACESIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SdispReplaceSigma = atof( value.c_str() );
			continue;
		}
		// if after secondary displacement replace mean
		key = "SDISPREPLACEMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SdispReplaceMean = atof( value.c_str() );
			continue;
		}
		// if after secondary displacement smooth sigma
		key = "SDISPLACESMOOTHSIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SdispSmoothSigma = atof( value.c_str() );
			continue;
		}
		// if after secondary displacement smooth mean
		key = "SDISPLACESMOOTHMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SdispSmoothMean = atof( value.c_str() );
			continue;
		}
		
		// if secondary strain error tollerance
		key = "SSTRAINERRORTOLLERANCE";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SstrainErrorToll = atof( value.c_str() );
			continue;
		}		
		// if after secondary strain replace sigma
		key = "SSTRAINREPLACESIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SstrainReplaceSigma = atof( value.c_str() );
			continue;
		}
		// if after secondary strain replace mean
		key = "SSTRAINREPLACEMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SstrainReplaceMean = atof( value.c_str() );
			continue;
		}
		// if after secondary strain smooth sigma
		key = "SSTRAINSMOOTHSIGMA";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SstrainSmoothSigma = atof( value.c_str() );
			continue;
		}
		// if after secondary strain smooth mean
		key = "SSTRAINSMOOTHMEAN";
		if ( !cLine.compare(0,key.size(),key) ){
			value.assign(cLine,key.size()+1,511);
			this->m_SstrainSmoothMean = atof( value.c_str() );
			continue;
		}
		
		// if no keys are found, make sure to give the loop a chance to exit
		if ( configFileInput.eof() ) {continue;}
		
		std::cout<<"Unknown key, value pair detected.  Please check the configuration file."<<std::endl;
		std::cout<<"Erronious line:"<<std::endl;
		std::cout<<cLine<<std::endl<<std::endl;
		return 1;
	}
	
	return 0;	
}

/** A function to read and set the fixed image file */
void ReadFixedImage()
{
	typedef itk::ImageFileReader<FixedImageType>		FixedImageReaderType;
	typename FixedImageReaderType::Pointer reader = FixedImageReaderType::New();
	reader->SetFileName(this->m_fixedFileName);
	
	try{
		reader->Update();
	}
	catch( itk::ExceptionObject &err ){
		std::cout<<"Error reading fixed image."<<std::endl<<"Message: "<<std::endl;
		std::cout<<err<<std::endl;
		std::exit(1);
	}
	
	this->SetFixedImage( reader->GetOutput() );
}


/** A function to read and set the moving image file */
void ReadMovingImage()
{
	typedef itk::ImageFileReader<MovingImageType>		MovingImageReaderType;
	typename MovingImageReaderType::Pointer reader = MovingImageReaderType::New();
	reader->SetFileName(this->m_movingFileName);
	
	try{
		reader->Update();
	}
	catch( itk::ExceptionObject &err ){
		std::cout<<"Error reading moving image."<<std::endl<<"Message: "<<std::endl;
		std::cout<<err<<std::endl;
		std::exit(1);
	}
	
	this->SetMovingImage( reader->GetOutput() );
}

void ReadMeshFile()
{
	if ( !this->m_meshFileName.compare(this->m_meshFileName.size()-3,3,"vtk") ){
		this->ReadVTKMesh( this->m_meshFileName );
		this->m_RestartFile = 1;				
	};
	if ( !this->m_meshFileName.compare(this->m_meshFileName.size()-3,3,"msh") ){this->ReadMeshFromGmshFile( this->m_meshFileName );}
}

/** A function to print the current setup of the DVC.  Returns a string
 * that can be written to the log file or to another output. */
std::string PrintConfiguration()
{
	std::stringstream outputText("");
	
	outputText<<"FIXEDIMAGEFILE="<<this->m_fixedFileName<<std::endl;
	outputText<<"MOVINGIMAGEFILE="<<this->m_movingFileName<<std::endl;
	outputText<<"MESHFILENAME="<<this->m_meshFileName<<std::endl;
	outputText<<"OUTPUTFOLDER="<<this->m_outputDirectory<<std::endl;
	outputText<<"IRRADIUS="<<this->GetInterrogationRegionRadius()<<std::endl;
	outputText<<"NTHREADS="<<this->GetRegistrationMethod()->GetNumberOfThreads()<<std::endl;
	outputText<<"GLOBALMAXSTEP="<<this->m_GlobalMaxStep<<std::endl;
	outputText<<"GLOBALMINSTEP="<<this->m_GlobalMinStep<<std::endl;
	outputText<<"INITIALDVCMAXSTEP="<<this->m_InitialDVCMaxStep<<std::endl;
	outputText<<"INITIALDVCMINSTEP="<<this->m_InitialDVCMinStep<<std::endl;
	outputText<<"SECONDARYDVCMAXSTEP="<<this->m_SecondaryDVCMaxStep<<std::endl;
	outputText<<"SECONDARYDVCMINSTEP="<<this->m_SecondaryDVCMinStep<<std::endl;
	//~ outputText<<"TERTIARYDVCMAXSTEP="<<this->m_TertiaryDVCMaxStep<<std::endl;
	//~ outputText<<"TERTIARYDVCMAXSTEP="<<this->m_TertiaryDVCMinStep<<std::endl;
	outputText<<"PERFORMSECONDARYDVC="<<this->m_SecondaryDVC<<std::endl;
	//~ outputText<<"PERFORMTERTIARTYDVC="<<this->m_TertiaryDVC<<std::endl;
	outputText<<"IDISPLACEMENTERRORTOLLERANCE="<<this->m_IdispErrorToll;
	outputText<<"IDISPREPLACESIGMA="<<this->m_IdispReplaceSigma<<std::endl;
	outputText<<"IDISPREPLACEMEAN="<<this->m_IdispReplaceMean<<std::endl;
	outputText<<"IDISPLACESMOOTHSIGMA="<<this->m_IdispSmoothSigma<<std::endl;
	outputText<<"IDISPLACESMOOTHMEAN="<<this->m_IdispSmoothMean<<std::endl;
	outputText<<"ISTRAINERRORTOLLERANCE="<<this->m_IstrainErrorToll;
	outputText<<"ISTRAINREPLACESIGMA="<<this->m_IstrainReplaceSigma<<std::endl;
	outputText<<"ISTRAINREPLACEMEAN="<<this->m_IstrainReplaceMean<<std::endl;
	outputText<<"ISTRAINSMOOTHSIGMA="<<this->m_IstrainSmoothSigma<<std::endl;
	outputText<<"ISTRAINSMOOTHMEAN="<<this->m_IstrainSmoothMean<<std::endl;
	outputText<<"SDISPLACEMENTERRORTOLLERANCE="<<this->m_SdispErrorToll;
	outputText<<"SDISPREPLACESIGMA="<<this->m_SdispReplaceSigma<<std::endl;
	outputText<<"SDISPREPLACEMEAN="<<this->m_SdispReplaceMean<<std::endl;
	outputText<<"SDISPLACESMOOTHSIGMA="<<this->m_SdispSmoothSigma<<std::endl;
	outputText<<"SDISPLACESMOOTHMEAN="<<this->m_SdispSmoothMean<<std::endl;
	outputText<<"SSTRAINERRORTOLLERANCE="<<this->m_SstrainErrorToll;
	outputText<<"SSTRAINREPLACESIGMA="<<this->m_SstrainReplaceSigma<<std::endl;
	outputText<<"SSTRAINREPLACEMEAN="<<this->m_SstrainReplaceMean<<std::endl;
	outputText<<"SSTRAINSMOOTHSIGMA="<<this->m_SstrainSmoothSigma<<std::endl;
	outputText<<"SSTRAINSMOOTHMEAN="<<this->m_SstrainSmoothMean<<std::endl;
	outputText<<std::endl;
	
	return outputText.str();
}

/** A function to get the current time */
std::string GetTime()
{
	struct tm * timeValue;
	std::time_t rawTime;
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	std::string currentTime = std::asctime( timeValue );
	
	return currentTime;
}

/** Setup the gobal registration.  This will use a linear interpolator.*/
void SetUpGlobalRegistration()
{
	/** get the registration method from the DVC algorithm */
	typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodPointer		registration = this->GetRegistrationMethod();
	typename DICMesh<FixedImageType,MovingImageType>::TransformTypePointer				transform = this->GetTransform();
	typename DICMesh<FixedImageType,MovingImageType>::OptimizerTypePointer				optimizer = this->GetOptimizer();
	typedef typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType	ParametersType;
	
	typedef itk::LinearInterpolateImageFunction< MovingImageType, double >			InterpolatorType;
	typename InterpolatorType::Pointer	linearInterpolator = InterpolatorType::New();
	registration->SetInterpolator( linearInterpolator );
	
	transform->SetIdentity();
	ParametersType	initialParameters = transform->GetParameters();
	registration->SetInitialTransformParameters( initialParameters );

	this->GetObserver()->SetLogfileName( this->GetLogfileName() );
	optimizer->AddObserver( itk::IterationEvent(), this->GetObserver() );
	
	/** The rotation part of the optimization is expected to be small 
	 * and it is more sensitive.  Use rotation values 3% of the translations. 
	 * There is some help on the ITK wiki for defining them:
	 * http://www.vtk.org/Wiki/ITK/ImageRegistration#Optimizers_2 */
	typedef typename DICMesh<FixedImageType,MovingImageType>::OptimizerType::ScalesType		OptimizerScalesType;
	OptimizerScalesType optScales( transform->GetNumberOfParameters() );
	optScales[0]	= 100;
	optScales[1]	= 100;
	optScales[2]	= 100;
	optScales[3]	= 100;
	optScales[4]	= 100;
	optScales[5]	= 100;
	optScales[6]	= 100;
	optScales[7]	= 100;
	optScales[8]	= 100;
	optScales[9]	= .05;
	optScales[10]	= .05;
	optScales[11]	= .05;
	optScales[12]	= .05;
	optScales[13]	= .05;
	optScales[14]	= .05;
	optimizer->SetScales( optScales );
	
	optimizer->SetMaximumStepLength( this->m_GlobalMaxStep );
	optimizer->SetMinimumStepLength( this->m_GlobalMinStep );
}

void SetupInitialDVCRegistration()
{
	/** get the registration method from the DVC algorithm */
	typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodPointer		registration = this->GetRegistrationMethod();
	typename DICMesh<FixedImageType,MovingImageType>::TransformTypePointer				transform = this->GetTransform();
	typename DICMesh<FixedImageType,MovingImageType>::OptimizerTypePointer				optimizer = this->GetOptimizer();
	typedef typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType	ParametersType;
	
	typedef itk::BSplineInterpolateImageFunction<FixedImageType, double, double>  BSplineInterpolatorType;
	typename BSplineInterpolatorType::Pointer	bSplineInterpolator = BSplineInterpolatorType::New();
	bSplineInterpolator->SetSplineOrder( 4 );
	registration->SetInterpolator( bSplineInterpolator );
	
	this->GetObserver()->SetLogfileName( this->GetLogfileName() );
	optimizer->AddObserver( itk::IterationEvent(), this->GetObserver() );
	
	/** The rotation part of the optimization is expected to be small 
	 * and it is more sensitive.  Use rotation values 3% of the translations. 
	 * There is some help on the ITK wiki for defining them:
	 * http://www.vtk.org/Wiki/ITK/ImageRegistration#Optimizers_2 */
	typedef typename DICMesh<FixedImageType,MovingImageType>::OptimizerType::ScalesType		OptimizerScalesType;
	OptimizerScalesType optScales( transform->GetNumberOfParameters() );
	optScales[0]	= 100;
	optScales[1]	= 100;
	optScales[2]	= 100;
	optScales[3]	= 100;
	optScales[4]	= 100;
	optScales[5]	= 100;
	optScales[6]	= 100;
	optScales[7]	= 100;
	optScales[8]	= 100;
	optScales[9]	= .05;
	optScales[10]	= .05;
	optScales[11]	= .05;
	optScales[12]	= .05;
	optScales[13]	= .05;
	optScales[14]	= .05;
	optimizer->SetScales( optScales );
	
	optimizer->SetMaximumStepLength( this->m_InitialDVCMaxStep );
	optimizer->SetMinimumStepLength( this->m_InitialDVCMinStep );
	
	this->CalculateInitialFixedImageRegionList();
	this->CalculateInitialMovingImageRegionList();
}

void SetupSecondaryDVCRegistration()
{
	/** get the registration method from the DVC algorithm */
	typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodPointer		registration = this->GetRegistrationMethod();
	typename DICMesh<FixedImageType,MovingImageType>::TransformTypePointer				transform = this->GetTransform();
	typename DICMesh<FixedImageType,MovingImageType>::OptimizerTypePointer				optimizer = this->GetOptimizer();
	typedef typename DICMesh<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType	ParametersType;
	
	typedef itk::BSplineInterpolateImageFunction<FixedImageType, double, double>  BSplineInterpolatorType;
	typename BSplineInterpolatorType::Pointer	bSplineInterpolator = BSplineInterpolatorType::New();
	bSplineInterpolator->SetSplineOrder( 4 );
	registration->SetInterpolator( bSplineInterpolator );
	
	this->GetObserver()->SetLogfileName( this->GetLogfileName() );
	optimizer->AddObserver( itk::IterationEvent(), this->GetObserver() );
	
	/** The rotation part of the optimization is expected to be small 
	 * and it is more sensitive.  Use rotation values 3% of the translations. 
	 * There is some help on the ITK wiki for defining them:
	 * http://www.vtk.org/Wiki/ITK/ImageRegistration#Optimizers_2 */
	typedef typename DICMesh<FixedImageType,MovingImageType>::OptimizerType::ScalesType		OptimizerScalesType;
	OptimizerScalesType optScales( transform->GetNumberOfParameters() );
	optScales[0]	= 100;
	optScales[1]	= 100;
	optScales[2]	= 100;
	optScales[3]	= 100;
	optScales[4]	= 100;
	optScales[5]	= 100;
	optScales[6]	= 100;
	optScales[7]	= 100;
	optScales[8]	= 100;
	optScales[9]	= .05;
	optScales[10]	= .05;
	optScales[11]	= .05;
	optScales[12]	= .05;
	optScales[13]	= .05;
	optScales[14]	= .05;
	optimizer->SetScales( optScales );
	
	optimizer->SetMaximumStepLength( this->m_SecondaryDVCMaxStep );
	optimizer->SetMinimumStepLength( this->m_SecondaryDVCMinStep );
	
	this->CalculateInitialFixedImageRegionList();
	this->CalculateInitialMovingImageRegionList();
}

bool RestartAnalysis()
{
	return m_RestartFile;
}

CommandIterationUpdate::Pointer GetObserver()
{
	return this->m_observer;
}

bool PerformSecondaryDVC()
{
	return this->m_SecondaryDVC;
}

//~ bool PerformTertiaryDVC()
//~ {
	//~ return this->m_TertiaryDVC;
//~ }

void ReplaceDisplacementBadPixelsAfterInitialDVC()
{
	this->SetDisplacementErrorTolerance( this->m_IdispErrorToll );
	
	vtkSmartPointer<vtkIdList> replacedPixels = vtkSmartPointer<vtkIdList>::New();
	this->ReplaceBadDisplacementPixels( this->m_IdispReplaceSigma, this->m_IdispReplaceMean, replacedPixels );
	std::stringstream msg("");
	for ( int i = 0; i < replacedPixels->GetNumberOfIds(); ++i){
		msg <<"Pixel "<<replacedPixels->GetId( i )<<" replaced."<<std::endl;
	}
	this->WriteToLogfile( msg.str() );
}

void SmoothDisplacementAfterInitialDVC()
{
	this->SetDisplacementErrorTolerance( this->m_IdispErrorToll );
	
	this->DisplacementWeightedMovingAverageFilter( this->m_IdispSmoothSigma, this->m_IdispSmoothMean );
}

void ReplaceStrainBadPixelsAfterInitialDVC()
{
	this->SetStrainErrorTolerance( this->m_IstrainErrorToll );
	
	vtkSmartPointer<vtkIdList> replacedPixels = vtkSmartPointer<vtkIdList>::New();
	this->ReplaceBadStrainPixels( this->m_IstrainReplaceSigma, this->m_IstrainReplaceMean, replacedPixels );
	std::stringstream msg("");
	for ( int i = 0; i < replacedPixels->GetNumberOfIds(); ++i){
		msg <<"Pixel "<<replacedPixels->GetId( i )<<" replaced."<<std::endl;
	}
	this->WriteToLogfile( msg.str() );
}

void SmoothStrainAfterInitialDVC()
{
	this->SetStrainErrorTolerance( this->m_IstrainErrorToll );
	
	this->StrainWeightedMovingAverageFilter( this->m_IstrainSmoothSigma, this->m_IstrainSmoothMean );
}

void ReplaceDisplacementBadPixelsAfterSecondDVC()
{
	this->SetDisplacementErrorTolerance( this->m_SdispErrorToll );

	vtkSmartPointer<vtkIdList> replacedPixels = vtkSmartPointer<vtkIdList>::New();
	this->ReplaceBadDisplacementPixels( this->m_SdispReplaceSigma, this->m_SdispReplaceMean, replacedPixels );
	std::stringstream msg("");
	for ( int i = 0; i < replacedPixels->GetNumberOfIds(); ++i){
		msg <<"Pixel "<<replacedPixels->GetId( i )<<" replaced."<<std::endl;
	}
	this->WriteToLogfile( msg.str() );
}

void SmoothDisplacementAfterSecondDVC()
{
	this->SetDisplacementErrorTolerance( this->m_SdispErrorToll );
	
	this->DisplacementWeightedMovingAverageFilter( this->m_SdispSmoothSigma, this->m_SdispSmoothMean );
}

void ReplaceStrainBadPixelsAfterSecondDVC()
{
	this->SetStrainErrorTolerance( this->m_SstrainErrorToll );
	
	vtkSmartPointer<vtkIdList> replacedPixels = vtkSmartPointer<vtkIdList>::New();
	this->ReplaceBadStrainPixels( this->m_SstrainReplaceSigma, this->m_SstrainReplaceMean, replacedPixels );
	std::stringstream msg("");
	for ( int i = 0; i < replacedPixels->GetNumberOfIds(); ++i){
		msg <<"Pixel "<<replacedPixels->GetId( i )<<" replaced."<<std::endl;
	}
	this->WriteToLogfile( msg.str() );
}

void SmoothStrainAfterSecondDVC()
{
	this->SetStrainErrorTolerance( this->m_SstrainErrorToll );
	
	this->StrainWeightedMovingAverageFilter( this->m_SstrainSmoothSigma, this->m_SstrainSmoothMean );
}

private:

ConfigurationFileNameType	m_configFileName;

// Global Registration Parameters
StepLengthType			m_GlobalMinStep;
StepLengthType			m_GlobalMaxStep;

// Initial DVC Parameters
StepLengthType			m_InitialDVCMinStep;
StepLengthType			m_InitialDVCMaxStep;

// Secondary DVC Parameters
StepLengthType			m_SecondaryDVCMinStep;
StepLengthType			m_SecondaryDVCMaxStep;

// Tertiary DVC Parameters
StepLengthType			m_TertiaryDVCMinStep;
StepLengthType			m_TertiaryDVCMaxStep;

// Error detection and correction after initial
double					m_IdispErrorToll;
// Displacement replacement
double					m_IdispReplaceSigma;
double					m_IdispReplaceMean;
// Displacement smoothing
double					m_IdispSmoothSigma;
double					m_IdispSmoothMean;

double					m_IstrainErrorToll;
// Strain repalcement
double					m_IstrainReplaceSigma;
double					m_IstrainReplaceMean;
// Strain smoothing
double					m_IstrainSmoothSigma;
double					m_IstrainSmoothMean;

// Error detection and correction after secondary
double_t				m_SdispErrorToll;
// Displacement replacement
double					m_SdispReplaceSigma;
double					m_SdispReplaceMean;
// Displacement smoothing
double					m_SdispSmoothSigma;
double					m_SdispSmoothMean;

double					m_SstrainErrorToll;
// Strain replacement
double					m_SstrainReplaceSigma;
double					m_SstrainReplaceMean;
// Strain smoothing
double					m_SstrainSmoothSigma;
double					m_SstrainSmoothMean;

// switch for number of DVCs
bool					m_SecondaryDVC;
bool					m_TertiaryDVC;

// restart file indicator
bool					m_RestartFile;

// image file names
std::string				m_fixedFileName;
std::string				m_movingFileName;
std::string				m_meshFileName;
std::string				m_outputDirectory;

// registration observer
CommandIterationUpdate::Pointer m_observer;

};

#endif //ANALYZEDVC_H
