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

#include "AnalyzeDVC.cxx"

int main(int argc, char **argv)
{
	if ( argc < 2 || argc > 2 )
	{
		std::cerr<<"Improper arguments!"<<std::endl;
		std::cerr<<"Usage:"<<std::endl;
		std::cerr<<argv[0]<<" ConfigureationFile"<<std::endl;
		return EXIT_FAILURE;
	}
	
	/** define the images types*/
	typedef	short			ImagePixelType;
	const unsigned int	dimension = 3;
	typedef itk::Image< ImagePixelType, dimension >		FixedImageType;
	typedef	itk::Image< ImagePixelType, dimension >		MovingImageType;
	
	/** Create the AnalyzeDVC fileter */
	typedef AnalyzeDVC< FixedImageType, MovingImageType >		DVCType;
	DVCType		*dvcMethod = new DVCType;
	
	/** input the config file */
	std::string configFile = argv[1];
	dvcMethod->SetConfigurationFile( configFile );
	
	/** process the config file */
	bool readFail = dvcMethod->ReadConfigureationFile();
	if ( readFail ){
		std::cout<<"There was an error in the configuration file."<<std::endl;
		std::cout<<"Aborting."<<std::endl<<std::endl;
		return 1;
	}
	
	std::string message;
	{
		std::string argv0 = argv[0];
		std::string argv1 = argv[1];
		message = argv0+" "+argv1+"\n";
	}
	dvcMethod->WriteToLogfile( message );
	
	dvcMethod->WriteToLogfile( dvcMethod->PrintConfiguration() );
	
	message = "Algorithm Started at: "+dvcMethod->GetTime();
	dvcMethod->WriteToLogfile( message );
	
	/** Read the input files */
	message = "Reading fixed image.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReadFixedImage();
	
	message = "Reading moving image.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReadMovingImage();
	
	message = "Reading mesh file.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReadMeshFile();
	
	if ( !dvcMethod->RestartAnalysis() ){
		// setup the global registration
		dvcMethod->SetUpGlobalRegistration();
		// perform global registraion
		dvcMethod->GlobalRegistration();
		message = "Global registration completed at: "+dvcMethod->GetTime();
		dvcMethod->WriteToLogfile( message );
		
		// setup the initial DVC
		dvcMethod->SetupInitialDVCRegistration();
		// perform the initial DVC
		dvcMethod->ExecuteDIC();
		message = "Initial DVC completed at: "+dvcMethod->GetTime();
		dvcMethod->WriteToLogfile( message );
	}
	
	if ( dvcMethod->RestartAnalysis() ){
		message = "Restarting by skipping the global registration and inital DIC.\n The input image will still be checke for errors and smoothed.\n";
		dvcMethod->WriteToLogfile( message );
	}
	
	message = "Calculating Strains.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->GetStrains();

	message = "Calculating Principal Strains.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->GetPrincipalStrains();
	
	message = "Writing initial DVC results image to "+dvcMethod->GetOutputDirectory()+"/InitialDVC.vtk";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->WriteMeshToVTKFile( dvcMethod->GetOutputDirectory()+"/InitialDVC.vtk" );
	
	message = "Removing and replacing bad displacemnet data points.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReplaceDisplacementBadPixelsAfterInitialDVC();
	
	message = "Smoothing displacement data.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->SmoothDisplacementAfterInitialDVC();
	
	message = "Removing and replacing bad strain data points.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReplaceStrainBadPixelsAfterInitialDVC();
	
	message = "Smoothing strain data.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->SmoothStrainAfterInitialDVC();	
	
	message = "Writing post processed initail DVC results image to "+dvcMethod->GetOutputDirectory()+"/PostProcessedInitialDVC.vtk";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->WriteMeshToVTKFile( dvcMethod->GetOutputDirectory()+"/PostProcessedInitialDVC.vtk" );
	
	// if no secondary DVC is requested, stop the analysis
	if ( !dvcMethod->PerformSecondaryDVC() ) {
		message = "Analysis completed at: "+dvcMethod->GetTime();
		dvcMethod->WriteToLogfile( message );
		return 0;
	}
	
	// setup the secondary DVC
	dvcMethod->SetupSecondaryDVCRegistration();
	message = "Starting second round DVC at: "+dvcMethod->GetTime();
	dvcMethod->WriteToLogfile( message );
	
	// second round DVC
	dvcMethod->ExecuteDIC();

	message = "Second round DVC completed at: "+dvcMethod->GetTime();
	dvcMethod->WriteToLogfile( message );
	
	message = "Calculating Strains.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->GetStrains();

	message = "Calculating Principal Strains.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->GetPrincipalStrains();
	
	message = "Writing second round DVC results image to "+dvcMethod->GetOutputDirectory()+"/SecondDVC.vtk";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->WriteMeshToVTKFile( dvcMethod->GetOutputDirectory()+"/SecondDVC.vtk" );
	
	message = "Removing and replacing bad displacemnet data points.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReplaceDisplacementBadPixelsAfterSecondDVC();
	
	message = "Smoothing displacement data.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->SmoothDisplacementAfterSecondDVC();
	
	message = "Removing and replacing bad strain data points.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->ReplaceStrainBadPixelsAfterSecondDVC();
	
	message = "Smoothing strain data.";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->SmoothStrainAfterSecondDVC();
	
	message = "Writing post processed initail DVC results image to "+dvcMethod->GetOutputDirectory()+"/PostProcessedSecondDVC.vtk";
	dvcMethod->WriteToLogfile( message );
	dvcMethod->WriteMeshToVTKFile( dvcMethod->GetOutputDirectory()+"/PostProcessedSecondDVC.vtk" );
	
	message = "Analysis completed at: "+dvcMethod->GetTime();
	dvcMethod->WriteToLogfile( message );
	
	return 0;
}
