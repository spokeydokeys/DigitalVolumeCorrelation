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


/** This file is used to test the functionality of new algorithms.  It
 * loads a generic vtkUnstructuredGrid into memeory and then will perform
 * operations frim the DICMesh class on it. */

int main(int argc, char** argv)
{
	if( argc < 4 || argc > 4 ){
		std::cout<<"Fatal Error: Incorrect Usage."<<std::endl;
		std::cout<<"Usage:"<<std::endl<<argv[0]<<" [Input Unstructured Grid] [Fixed Image] [Moving Image]"<<std::endl;
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

	
	FixedReaderType::Pointer fixedReader = FixedReaderType::New();
	MovingReaderType::Pointer movingReader = MovingReaderType::New();	
	
	fixedReader->SetFileName( argv[2] );
	movingReader->SetFileName( argv[3] );
	
	fixedReader->Update();
	movingReader->Update();
	
	// define and instantate the DICmethod
	typedef	DICMesh<FixedImageType, MovingImageType>	DICType;
	DICType *DICMethod = new DICType;
	
	/*std::string logfileName = "/home/seth/logfile.txt";
	DICMethod->SetLogfileName( logfileName );*/
	DICMethod->SetFixedImage(fixedReader->GetOutput() );
	DICMethod->SetMovingImage(movingReader->GetOutput() );
	
	DICMethod->SetInterrogationRegionRadius(14);
	/*std::string outputDir = "/home/seth";
	DICMethod->SetOuputDirectory( outputDir );*/
		
	// import a vtkUnstructuredGrid.
	vtkSmartPointer<vtkUnstructuredGridReader>		vtkReader= vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkReader->SetFileName( argv[1] );
	vtkReader->Update();
	DICMethod->SetDataImage( vtkReader->GetOutput() );
	
	// Test the algorithm here.
	DICMethod->CreateNewRegionListFromBadPixels();
	//DICMethod->ExecuteDIC();
	//DICMethod->GetStrains();
	//DICMethod->GetPrincipalStrains();	
	
	// Write the results of the test
//	std::string outFile = argv[2];
//	DICMethod->WriteMeshToVTKFile( outFile );
	
	return 0;
}
