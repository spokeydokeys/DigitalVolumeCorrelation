//      Test.cxx
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

int main(int argc, char **argv)
{
	if ( argc < 6 || argc > 6 )
	{
		std::cerr<<"Improper arguments!"<<std::endl;
		std::cerr<<"Usage:"<<std::endl;
		std::cerr<<argv[0]<<" FixedImage MovingImage GmshFile OutputDirectory IRRegionRadius"<<std::endl;
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
	
	/* This block would be used to import a vtk file, rather than read a gmsh file. */
	//~ vtkSmartPointer<vtkUnstructuredGridReader>		vtkReader= vtkSmartPointer<vtkUnstructuredGridReader>::New();
	//~ vtkReader->SetFileName( argv[3] );
	//~ vtkReader->Update();
	//~ DICMethod->SetDataImage( vtkReader->GetOutput() );
	
	
	DICMethod->ReadMeshFromGmshFile( argv[3] );
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
	
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();	// thes lines will make the optimizer print out its
	observer->SetLogfileName( DICMethod->GetLogfileName() );					// displacement as it goes.  They can be removed.
	optimizer->AddObserver( itk::IterationEvent(), observer );
	
	// The rotation part of the optimization is expected to be small 
	// and it is more sensitive.  Use rotation values 3% of the translations.
	typedef DICType::OptimizerType::ScalesType		OptimizerScalesType;
	OptimizerScalesType optScales( transform->GetNumberOfParameters() );	// optimizer scales must be chosen carefully.  These were based off of
	optScales[0] = 10;														// prelininary tests on my data.  There is some help on the ITK wiki
	optScales[1] = 10;														// for defining them:
	optScales[2] = 10;														// http://www.vtk.org/Wiki/ITK/ImageRegistration#Optimizers_2
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
	optimizer->SetMaximumStepLength(.5); // large steps for the global registration (based on visual alignment in ParaView)
	optimizer->SetMinimumStepLength(0.041); // low tolerance for the global registration

	{ // scoping the down-sampling will make the Shrink filters be automatically distroyed (freeing memory)
		typedef itk::ShrinkImageFilter< FixedImageType, FixedImageType > ResamplerType;
		ResamplerType::Pointer	fixedResampler = ResamplerType::New();
		ResamplerType::Pointer	movingResampler = ResamplerType::New();		
		fixedResampler->SetInput( fixedReader->GetOutput() );
		movingResampler->SetInput( movingReader->GetOutput() );
		fixedResampler->SetShrinkFactors(3);
		movingResampler->SetShrinkFactors(3);
		fixedResampler->SetNumberOfThreads(8);
		movingResampler->SetNumberOfThreads(8);
		msg.str("");
		msg <<"Resampling for global registration"<<std::endl<<std::endl;
		DICMethod->WriteToLogfile(msg.str());
		fixedResampler->Update();
		movingResampler->Update();
		// global registration - rotation is centred on the body
		registration->SetFixedImage( fixedResampler->GetOutput() );
		registration->SetMovingImage( movingResampler->GetOutput() );
		MovingImageType::PointType	movingOrigin = movingResampler->GetOutput()->GetOrigin();
		MovingImageType::SpacingType movingSpacing = movingResampler->GetOutput()->GetSpacing();
		MovingImageType::SizeType	movingSize = movingResampler->GetOutput()->GetLargestPossibleRegion().GetSize();
		MovingImageType::PointType	centerPt;
		centerPt[0] = movingOrigin[0] + movingSpacing[0]*(movingSize[0] >> 2);
		centerPt[1] = movingOrigin[1] + movingSpacing[1]*(movingSize[1] >> 2);
		centerPt[2] = movingOrigin[2] + movingSpacing[2]*(movingSize[2] >> 2);
		ParametersType initalParams = registration->GetInitialTransformParameters();
		initialParameters[3] = centerPt[0];
		initialParameters[4] = centerPt[1];
		initialParameters[5] = centerPt[2];
		registration->SetInitialTransformParameters( initialParameters );
		// restrict the Global registration to just the bounding box of the mesh
		double *meshBBox = new double[6];
		meshBBox = DICMethod->GetDataImage()->GetBounds();
		FixedImageType::PointType meshMinPt;
		meshMinPt[0] = *meshBBox;
		meshMinPt[1] = *(meshBBox+2);
		meshMinPt[2] = *(meshBBox+4);
		double meshSize[3];
		meshSize[0] = *(meshBBox+1)-*meshBBox;
		meshSize[1] = *(meshBBox+3)-*(meshBBox+2);
		meshSize[2] = *(meshBBox+5)-*(meshBBox+4);
		FixedImageType::IndexType fixedImageROIStart;
		fixedResampler->GetOutput()->TransformPhysicalPointToIndex(meshMinPt,fixedImageROIStart);
		FixedImageType::SpacingType fixedSpacing = fixedResampler->GetOutput()->GetSpacing();
		FixedImageType::SizeType fixedImageROILengths;
		fixedImageROILengths[0] = (int)std::floor(meshSize[0]/fixedSpacing[0]);
		fixedImageROILengths[1] = (int)std::floor(meshSize[1]/fixedSpacing[1]);
		fixedImageROILengths[2] = (int)std::floor(meshSize[2]/fixedSpacing[2]);
		FixedImageType::RegionType fixedAnalysisRegion;
		fixedAnalysisRegion.SetIndex( fixedImageROIStart );
		fixedAnalysisRegion.SetSize( fixedImageROILengths );
		registration->SetFixedImageRegion( fixedAnalysisRegion );
		registration->SetFixedImageRegionDefined( true );
		msg.str("");
		msg << "Global registration in progress"<<std::endl<<std::endl;
		DICMethod->WriteToLogfile( msg.str() );
		registration->Update();
		registration->SetFixedImageRegionDefined( false );
	}
	
	msg.str("");
	msg << "Global Registration complete."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	double *globalRegResults = new double[3];
	ParametersType	finalParameters = registration->GetLastTransformParameters();
	msg.str("");
	msg << "Final Params:"<< finalParameters<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	*globalRegResults = finalParameters[6];
	*(globalRegResults + 1) = finalParameters[7];
	*(globalRegResults + 2) = finalParameters[8];	
	
	msg.str("");
	msg << "Global registration finished.\n Resulting displacement: ("<<*globalRegResults<<
		", "<<*(globalRegResults+1)<<", "<<*(globalRegResults+2)<<")"<<std::endl<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	optimizer->SetMaximumStepLength( 0.041 ); // smaller steps for the DIC
	optimizer->SetMinimumStepLength( 0.0005); // increased tolerace for the DIC.
	DICMethod->SetMeshInitialValue( globalRegResults );
	DICMethod->SetRegistrationMethod( registration );
	
	//~ {
		//~ 
		//~ DICMethod->CalculateInitialFixedImageRegionList();
		//~ DICMethod->CalculateInitialMovingImageRegionList();
		//~ 
		//~ DICType::MovingImageRegionListType *movingRegionList = DICMethod->GetMovingImageRegionList();
		//~ 
		//~ for( unsigned int k = 0; k < movingRegionList->size(); ++k ){
			//~ std::cout<<"Region "<<k<<std::endl;
			//~ DICType::MovingImageType::RegionType *movingRegion = DICMethod->GetMovingImageRegionFromIndex( k );
			//~ DICType::FixedImageType::RegionType *fixedRegion = DICMethod->GetFixedImageRegionFromIndex( k );
			//~ std::cout<<"Moving Region: "<<*movingRegion<<std::endl;
			//~ std::cout<<"Fixed Region: "<<*fixedRegion<<std::endl;
		//~ }
		//~ 
	//~ }	
	
	DICMethod->ExecuteDIC();
	msg.str("");
	msg << "Intitial DVC finished."<<std::endl<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	msg.str("");
	msg << "Checking for bad pixels."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->CreateNewRegionListFromBadPixels();
	
	msg.str("");
	msg << "Starting second round DVC."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->ExecuteDIC();
	
	msg.str("");
	msg << "Second round DVC complete."<<std::endl<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	msg.str("");
	msg << "Checking for bad pixels."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->CreateNewRegionListFromBadPixels();
	
	msg.str("");
	msg << "Starting third round DVC."<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	DICMethod->ExecuteDIC();
	
	msg.str("");
	msg << "Third round DVC complete."<<std::endl<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
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
	std::string outputFile = outputDir + "/DIC_Result.vtk";
	DICMethod->WriteMeshToVTKFile( outputFile );
	
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	msg.str("");
	msg << "Finishing DIC at: "<<std::asctime( timeValue )<<std::endl;
	DICMethod->WriteToLogfile( msg.str() );
	
	return 0;
}

