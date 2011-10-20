//      AnalyzeImages.cxx
//      
//      Copyright 2011 Seth Gilchrist <seth@mech.ubc.ca>
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


#include "itkShrinkImageFilter.h"
#include "DICMesh.cxx"
#include "itkImageFileWriter.h"
#include <vtkUnstructuredGridReader.h>
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLBFGSBOptimizer.h"
#include <itkBSplineInterpolateImageFunction.h>


// the following provides updates for the RegularStep optimizer
class RSOCommandIterationUpdate : public itk::Command
{
public:
typedef  RSOCommandIterationUpdate   Self;
typedef  itk::Command             Superclass;
typedef itk::SmartPointer<Self>  Pointer;
itkNewMacro( Self );

protected:
RSOCommandIterationUpdate() {};

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
	msg << "Iteration: " << optimizer->GetCurrentIteration() << "   ";
    msg << "Metric Value: " << optimizer->GetCachedValue() << "   ";
    msg << "Position: " << optimizer->GetCachedCurrentPosition() << "   ";
    msg << "Gradient Magnitude: " << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
	this->WriteToLogfile( msg.str() );
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



int main(int argc, char **argv)
{
	if ( argc < 6 || argc > 6 )
	{
		std::cerr<<"Improper arguments!"<<std::endl;
		std::cerr<<"Usage:"<<std::endl;
		std::cerr<<argv[0]<<" FixedImage MovingImage GmshFile OutputDirectory IRRegionRadius"<<std::endl<<std::endl;
		std::cerr<<"If a .vtk file is provided in place of a GmshFile it will be assumed that the analysis is a restart and the global and initial DVC will be skipped."<<std::endl;
		std::cerr<<"Aborting"<<std::endl;
		return EXIT_FAILURE;
	}
	
	std::stringstream msg("");
	
	// create the images
	typedef	short			ImagePixelType;
	const unsigned int	dimension = 3;
	typedef itk::Image< ImagePixelType, dimension >		FixedImageType;
	typedef	itk::Image< ImagePixelType, dimension >		MovingImageType;
	
	// creating a DICMesh filter
	typedef	DICMesh< FixedImageType, MovingImageType >		DICType;
	DICType 	*DICMethod = new DICType;
	std::string outputDir = argv[4];
	std::string logfile = outputDir+"/logfile.txt";
	DICMethod->SetLogfileName( logfile );
	DICMethod->SetOuputDirectory( outputDir );
	
	msg << argv[0] <<" '"<< argv[1]<<"' '"<<argv[2]<<"' '"<<argv[3]<<"' '"<<argv[4]<<"' "<<argv[5]<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );

	struct tm * timeValue;
	std::time_t rawTime;
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	msg.str("");
	msg << "Algorithm started at: "<<std::asctime( timeValue )<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	// create the image readers
	typedef itk::ImageFileReader< FixedImageType >		FixedImageReaderType;
	typedef	itk::ImageFileReader< MovingImageType >		MovingImageReaderType;
	FixedImageReaderType::Pointer	fixedReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer	movingReader = MovingImageReaderType::New();
	fixedReader->SetFileName( argv[1] );
	movingReader->SetFileName( argv[2] );
	msg.str("");
	msg << "Reading the images."<<std::endl<<std::endl;
	DICMethod->WriteToLogfile(msg.str());
	
	fixedReader->Update();
	movingReader->Update();
	
	DICMethod->SetFixedImage(fixedReader->GetOutput());
	DICMethod->SetMovingImage( movingReader->GetOutput() );
	
	std::string inputMesh = argv[3];
	bool restartFile = false;
	inputMesh.find("vtk",inputMesh.length()-5) == std::string::npos ? restartFile = false : restartFile = true;
	std::cout<<inputMesh<<"  "<<restartFile<<std::endl;
	
	if ( restartFile ){
		restartFile = true;
		msg.str("");
		msg<<"Using "<<inputMesh<<" as a restart file."<<std::endl<<std::endl;
		DICMethod->WriteToLogfile(msg.str());
	}
	
	// If the file is a restart file, read it in and set it as the data image
	if ( restartFile ){
		vtkSmartPointer<vtkUnstructuredGridReader>		vtkReader= vtkSmartPointer<vtkUnstructuredGridReader>::New();
		vtkReader->SetFileName( argv[3] );
		vtkReader->Update();
		DICMethod->SetDataImage( vtkReader->GetOutput() );
	}
	else{		
		DICMethod->ReadMeshFromGmshFile( argv[3] );
	}
	
	msg.str("");
	msg << "Images in memory."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	// Set the interrogation region radius, the region will be 2*radius+1 larger than this number
	DICMethod->SetInterrogationRegionRadius( atoi(argv[5]) );
	msg.str("");
	msg <<"Interrogation region radius set to: "<<DICMethod->GetInterrogationRegionRadius()<<std::endl<<std::endl;// confirm region radius
	DICMethod->WriteToLogfile(msg.str());
	
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
	
	/*NOCommandIterationUpdate::Pointer nOObserver = NOCommandIterationUpdate::New();	// thes lines will make the optimizer print out itsO
	nOObserver->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	optimizer->AddObserver( itk::IterationEvent(), nOObserver );*/
	
	RSOCommandIterationUpdate::Pointer observer = RSOCommandIterationUpdate::New();	// thes lines will make the optimizer print out its
	observer->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	optimizer->AddObserver( itk::IterationEvent(), observer );
	
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
	optimizer->SetScales( optScales );	
	
	// global registration is in-exact and only gives an estimate for 
	// the rest of the DIC.  Using a downsampled image increases the 
	// radius of convergence and speeds things up.
	optimizer->SetMaximumStepLength(0.010); // large steps for the global registration (based on visual alignment in ParaView)
	optimizer->SetMinimumStepLength(0.005); // low tolerance for the global registration*/
	if ( !restartFile ) DICMethod->GlobalRegistration();
	
	// speed things us for the actual registration
	optimizer->SetMaximumStepLength( 0.001 ); // smaller steps for the DIC
	optimizer->SetMinimumStepLength( 0.0005); // increased tolerace for the DIC.
	
	if ( !restartFile ){
		DICMethod->ExecuteDIC();
		msg.str("");
		msg << "Intitial DVC finished."<<std::endl<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
	
		msg.str("");
		msg << "Calculating Strains"<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
		DICMethod->GetStrains();

		msg.str("");
		msg <<"Calculating Principal Strains."<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
		DICMethod->GetPrincipalStrains();
		
		std::string debugFile = DICMethod->GetOutputDirectory() + "/AfterInitialDVC.vtk";
		DICMethod->WriteMeshToVTKFile( debugFile );
		

	
		msg.str("");
		msg <<"Starting with the second round DVC"<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
	}
	
	if ( restartFile ){
		msg.str("");
		msg<<"Restarting by skipping the global registration and inital DIC."<<std::endl;
		msg<<"The input image will still be smoothed before starting the second round DVC."<<std::endl<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
	}
	
	msg.str("");
	msg << "Smoothing the image." <<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->WeightedMovingAverageFilter(10, 2, 0);
	
	DICMethod->CalculateInitialFixedImageRegionList();
	DICMethod->CalculateInitialMovingImageRegionList();
	
	// second round DVC gets a newtonian style optimizer. These are more accurate, but have a smaller
	// radius of convergence, that is why it is not used for the inital DVC
	typedef itk::LBFGSBOptimizer NewtonOptimizerType;
	NewtonOptimizerType::Pointer newtonOptimizer = NewtonOptimizerType::New();
	
	NOCommandIterationUpdate::Pointer nOObserver = NOCommandIterationUpdate::New();	// thes lines will make the optimizer print out its
	nOObserver->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	newtonOptimizer->AddObserver( itk::IterationEvent(), nOObserver );	

	newtonOptimizer->SetCostFunctionConvergenceFactor( 1e12 );
	newtonOptimizer->SetProjectedGradientTolerance( 1e-6 );
	
	registration->SetOptimizer( newtonOptimizer );
	
	// second round DVC gets a bspline interpolator. These are more accurate, but have a huge computational
	// cost, so the linear is used to get close in the first round
	
	typedef itk::BSplineInterpolateImageFunction<FixedImageType, double, double>  BSplineInterpolatorType;
	BSplineInterpolatorType::Pointer	bSplineInterpolator = BSplineInterpolatorType::New();
	bSplineInterpolator->SetSplineOrder( 4 );
	
	registration->SetInterpolator( bSplineInterpolator );
	
	// Execute the second round DVC
	DICMethod->ExecuteDIC();
	//~ 
	msg.str("");
	msg<<"Second Round Finished"<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	//~ msg.str("");
	//~ msg << "Checking for bad pixels."<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->CreateNewRegionListFromBadPixels();
	//~ 
	//~ msg.str("");
	//~ msg << "Starting second round DVC."<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->ExecuteDIC();
	//~ 
	msg.str("");
	msg << "Second round DVC complete."<<std::endl<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	//~ 
	//~ msg << "Calculating Strains"<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->GetStrains();
//~ 
	//~ msg.str("");
	//~ msg <<"Calculating Principal Strains."<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->GetPrincipalStrains();
	//~ 
	//~ debugFile = DICMethod->GetOutputDirectory() + "/AfterSecondDVC.vtk";
	//~ DICMethod->WriteMeshToVTKFile( debugFile );
	//~ 
	//~ msg.str("");
	//~ msg << "Checking for bad pixels."<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->CreateNewRegionListFromBadPixels();
	//~ //optimizer->SetMaximumStepLength( 0.0005 );
	//~ //optimizer->SetMinimumStepLength( 0.00005 );
	//~ 
	//~ msg.str("");
	//~ msg << "Starting third round DVC."<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	//~ DICMethod->ExecuteDIC();
	//~ 
	//~ msg.str("");
	//~ msg << "Third round DVC complete."<<std::endl<<std::endl;
	//~ DICMethod->WriteToLogfile( msg.str() );
	
	msg.str("");
	msg << "Calculating Strains"<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->GetStrains();

	msg.str("");
	msg <<"Calculating Principal Strains."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->GetPrincipalStrains();
	
	msg.str("");
	msg << "Writing output file"<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	std::string outputFile = outputDir + "/Final_DIC_Result.vtk";
	DICMethod->WriteMeshToVTKFile( outputFile );
	
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	msg.str("");
	msg << "Finishing DIC at: "<<std::asctime( timeValue )<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	return 0;
}

