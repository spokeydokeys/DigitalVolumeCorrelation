//      C++ Source.cpp
//      
//      Copyright 2010 Seth Gilchrist <seth@mech.ubc.ca>
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


#include <iostream>
#include "DICMesh.cxx"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include <vtkUnstructuredGridReader.h>
#include "itkImage.h"
#include "itkLBFGSBOptimizer.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"


/** This file is used to test the functionality of new algorithms.  It
 * loads a generic vtkUnstructuredGrid into memeory and then will perform
 * operations frim the DICMesh class on it. */

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
	msg << optimizer->GetCurrentIteration() << " = "  << optimizer->GetValue() << " : " << optimizer->GetCurrentPosition() << std::endl;
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

// the following provides updates for the Newtonian optimizer
class NOCommandIterationUpdate : public itk::Command
{
public:
typedef  NOCommandIterationUpdate   Self;
typedef  itk::Command             Superclass;
typedef itk::SmartPointer<Self>  Pointer;
itkNewMacro( Self );

protected:
NOCommandIterationUpdate() {};

public:

typedef itk::LBFGSBOptimizer     OptimizerType;
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
	//~ msg << optimizer->GetCurrentIteration() << " = "  << optimizer->GetValue() << " : " << optimizer->GetCurrentPosition() << std::endl;
	msg << optimizer->GetCurrentIteration() << "   ";
    msg << optimizer->GetCachedValue() << "   ";
    msg << optimizer->GetCachedCurrentPosition() << "   ";
    msg << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
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

int main(int argc, char** argv)
{
	if( argc < 5 || argc > 5 ){
		std::cout<<"Fatal Error: Incorrect Usage."<<std::endl;
		std::cout<<"Usage:"<<std::endl<<argv[0]<<" [Input Unstructured Grid] [Fixed Image] [Moving Image] [Output Path]"<<std::endl;
		std::cout<<"Exiting"<<std::endl;
		return EXIT_FAILURE;
	}
	
	// define image types.  needed to define the DICMethod.
	typedef	short		ImagePixelType;
	const unsigned int	dimension = 3;
	typedef itk::Image< ImagePixelType, dimension >		FixedImageType;
	typedef	itk::Image< ImagePixelType, dimension >		MovingImageType;
	typedef itk::ImageFileReader< FixedImageType >		FixedReaderType;
	typedef itk::ImageFileReader< MovingImageType >		MovingReaderType;
	
	// define and instantate the DICmethod
	typedef	DICMesh<FixedImageType, MovingImageType>	DICType;
	DICType *DICMethod = new DICType;
	std::string outputDir = argv[4];
	std::string logfile = outputDir+"/logfile.txt";
	DICMethod->SetLogfileName( logfile );
	DICMethod->SetOuputDirectory( outputDir );
	
	
	FixedReaderType::Pointer fixedReader = FixedReaderType::New();
	MovingReaderType::Pointer movingReader = MovingReaderType::New();	
	fixedReader->SetFileName( argv[2] );
	movingReader->SetFileName( argv[3] );
	fixedReader->Update();
	movingReader->Update();
	
	/*std::string logfileName = "/home/seth/logfile.txt";*/
	DICMethod->SetFixedImage(fixedReader->GetOutput() );
	DICMethod->SetMovingImage(movingReader->GetOutput() );
	
	
	DICMethod->SetInterrogationRegionRadius(36);
		
	// import a vtkUnstructuredGrid.
	vtkSmartPointer<vtkUnstructuredGridReader>		vtkReader= vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkReader->SetFileName( argv[1] );
	vtkReader->Update();
	DICMethod->SetDataImage( vtkReader->GetOutput() );
	
	// get the registration method from the DVC algorithm
	DICType::ImageRegistrationMethodPointer		registration = DICMethod->GetRegistrationMethod();
	DICType::TransformTypePointer				transform = DICMethod->GetTransform();
	DICType::OptimizerTypePointer				optimizer = DICMethod->GetOptimizer();
	typedef DICType::ImageRegistrationMethodType::ParametersType	ParametersType;
	
	// Setup the registration
	registration->SetNumberOfThreads(	8	);
	transform->SetIdentity();
	ParametersType	initialParameters = transform->GetParameters();
	registration->SetInitialTransformParameters( initialParameters );
	
	//~ CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();	// thes lines will make the optimizer print out its
	//~ observer->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	//~ optimizer->AddObserver( itk::IterationEvent(), observer );
	
	// The rotation part of the optimization is expected to be small 
	// and it is more sensitive.  Use rotation values 3% of the translations.
	typedef DICType::OptimizerType::ScalesType		OptimizerScalesType;
	OptimizerScalesType optScales( transform->GetNumberOfParameters() );	// optimizer scales must be chosen carefully.  These were based off of
	optScales[0] = 100;														// prelininary tests on my data.  There is some help on the ITK wiki
	optScales[1] = 100;														// for defining them:
	optScales[2] = 100;														// http://www.vtk.org/Wiki/ITK/ImageRegistration#Optimizers_2
	optScales[3] = .05;
	optScales[4] = .05;
	optScales[5] = .05;
	optScales[6] = .05;
	optScales[7] = .05;
	optScales[8] = .05;
	//~ optimizer->SetScales( optScales );	
	
	//~ // speed things us for the actual registration
	//~ optimizer->SetMaximumStepLength( 0.041 ); // smaller steps for the DIC
	//~ optimizer->SetMinimumStepLength( 0.0005); // increased tolerace for the DIC.
	
	typedef itk::LBFGSBOptimizer NewtonOptimizerType;
	NewtonOptimizerType::Pointer newtonOptimizer = NewtonOptimizerType::New();
	//~ newtonOptimizer->SetScales( optScales );
	//~ newtonOptimizer->SetGradientConvergenceTolerance( 1e-10 );
	//~ newtonOptimizer->SetDefaultStepLength( .00041 );
	//newtonOptimizer->SetCostFunction( registration->GetMetric() );
	
	typedef itk::NormalizedCorrelationImageToImageMetric<FixedImageType, MovingImageType> NormalizedMetricType;
	NormalizedMetricType::Pointer normalizedMetric = NormalizedMetricType::New();

	NOCommandIterationUpdate::Pointer nOObserver = NOCommandIterationUpdate::New();	// thes lines will make the optimizer print out itsO
	nOObserver->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	newtonOptimizer->AddObserver( itk::IterationEvent(), nOObserver );
	
	NewtonOptimizerType::BoundSelectionType boundSelect( transform->GetNumberOfParameters() );
	NewtonOptimizerType::BoundValueType upperBound( transform->GetNumberOfParameters() );
	NewtonOptimizerType::BoundValueType lowerBound( transform->GetNumberOfParameters() );

	boundSelect.Fill( 0 );
	upperBound.Fill( 0.0 );
	lowerBound.Fill( 0.0 );

	newtonOptimizer->SetBoundSelection( boundSelect );
	newtonOptimizer->SetUpperBound( upperBound );
	newtonOptimizer->SetLowerBound( lowerBound );

	newtonOptimizer->SetCostFunctionConvergenceFactor( 1.e7 );
	newtonOptimizer->SetProjectedGradientTolerance( 1e-6 );
	//~ newtonOptimizer->SetMaximumNumberOfIterations( 200 );
	//~ newtonOptimizer->SetMaximumNumberOfEvaluations( 30 );
	//~ newtonOptimizer->SetMaximumNumberOfCorrections( 5 );
	
	//~ registration->SetMetric( normalizedMetric );
	registration->SetOptimizer( newtonOptimizer );
	
	DICMethod->CalculateInitialFixedImageRegionList();
	DICMethod->CalculateInitialMovingImageRegionList();
	
	DICMethod->ExecuteDIC();
	
	std::string outFile = outputDir + "/result2.vtk";	
	DICMethod->WriteMeshToVTKFile( outFile );	
	
	//~ DICMethod->CreateNewRegionListFromBadPixels();
	//DICMethod->ExecuteDIC();
	//DICMethod->GetStrains();
	//DICMethod->GetPrincipalStrains();	
	
	// Write the results of the test

	
	return 0;
}
