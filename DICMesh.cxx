//      DICMesh.cxx
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

#ifndef DICMESH_H
#define DICMESH_H

#include <cstring>
#include <ctime>
#include "DIC.cxx"
#include "itkMesh.h"
#include "itkTetrahedronCell.h"
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkPKdTree.h>
#include <vtkIdList.h>
#include <vtkCellDerivatives.h>
#include <vtkCellDataToPointData.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageData.h>

#include <itkImageFileWriter.h>
#include <itkShrinkImageFilter.h>

template<typename TFixedImage, typename TMovingImage>
class DICMesh : public DIC<TFixedImage, TMovingImage>
{
public:

/** Inherited typedefs from DIC.*/
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageType				FixedImageType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageType				MovingImageType;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImagePointer			FixedImagePointer;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImagePointer			MovingImagePointer;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageType::IndexType 	FixedImageIndexType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageType::IndexType	MovingImageIndexType;
typedef	typename DIC<TFixedImage,TMovingImage>::FixedImageRegionType		FixedImageRegionType;
typedef typename DIC<TFixedImage,TMovingImage>::MovingImageRegionType		MovingImageRegionType;
typedef typename DIC<TFixedImage,TMovingImage>::FixedImageRegionListType	FixedImageRegionListType;
typedef	typename DIC<TFixedImage,TMovingImage>::MovingImageRegionListType	MovingImageRegionListType;

typedef typename DIC<TFixedImage,TMovingImage>::ImageRegistrationMethodType	ImageRegistrationMethodType;
typedef typename ImageRegistrationMethodType::ParametersType				RegistrationParametersType;	
typedef typename DIC<TFixedImage,TMovingImage>::TransformInitializerType	TransformInitializerType;
typedef typename DIC<TFixedImage,TMovingImage>::TransformType				TransformType;

typedef		vtkSmartPointer<vtkUnstructuredGrid>	DataImagePointer;
typedef		vtkSmartPointer<vtkPoints>				DataImagePointsPointer;
typedef		vtkSmartPointer<vtkDoubleArray>			DataImagePixelPointer;



/* Methods. **/

/** Constructor **/
DICMesh()
{
	m_DataImage = 0; // must be provided by user or read in using the ReadMeshFromGmshFile method
	m_KDTree = vtkSmartPointer<vtkPKdTree>::New();
	m_errorRadius = 4; // initialize the error search radius to 3 image units
	m_displacementErrorTolerance = 2; // difference of a pixel from its neighbours to be considered erronious, in standard deviations from the mean
	m_strainErrorTolerance = 1;	
	m_pointsList = vtkSmartPointer<vtkIdList>::New(); // the points list for analysis
	m_maxMeticValue = -0.00; // TODO: make this setable using a method
	m_GlobalRegDownsampleValue = 3; // This value is the default downsample when preforming the global registration.
}

/** Destructor **/
~DICMesh() {}

/** This function will compile a list of valid moving image regions.*/
void CalculateInitialMovingImageRegionList()
{
	vtkIdType	numberOfNodes = m_DataImage->GetNumberOfPoints();
	this->m_MovingImageRegionList.clear();
	
	// visit every point in the image
	for ( int i = 0; i < numberOfNodes; ++i){ 
		double *currentLocation = new double[3];
		currentLocation = this->GetMovingImageRegionLocationFromIndex( i ); // get the current node location
		
		MovingImageRegionType *currentRegion = new MovingImageRegionType;
		this->GetMovingImageRegionFromLocation( currentRegion, currentLocation ); // get the region
		this->PushRegionOntoMovingImageRegionList( currentRegion ); // add the region to the moving image region list.
	}
}

/** Given an intex, this will calcualted the center of a valid moving image region.*/
double *GetMovingImageRegionLocationFromIndex( vtkIdType i)
{
	double *currentPoint	= new double[3];
	double *currentValue	= new double[3];
	double *currentLocation	= new double[3];
	m_DataImage->GetPoint( i, currentPoint ); // get the current point location
	this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( i, currentValue );
	*currentLocation		= *currentPoint		+ *currentValue; // set the current point to the original local to the disp
	*(currentLocation +1)	= *(currentPoint +1)+ *(currentValue +1);
	*(currentLocation +2)	= *(currentPoint +2)+ *(currentValue +2);
	
	return currentLocation;
}

/** A function to calculate the initial fixed image region list. */
void CalculateInitialFixedImageRegionList()
{
	this->m_FixedImageRegionList.clear();
	this->m_pointsList->Reset();
	
	for (int i = 0; i < this->m_DataImage->GetNumberOfPoints(); ++i){
		double *currentLocation = new double[3];
		m_DataImage->GetPoint( i, currentLocation ); // the fixed image region is always centered on the current mesh point
	
		FixedImageRegionType	*currentRegion = new FixedImageRegionType;
		this->GetFixedImageRegionFromLocation( currentRegion, currentLocation );
		this->PushRegionOntoFixedImageRegionList( currentRegion );
		this->m_pointsList->InsertNextId( i );
	}
}

/** A function to read a gmsh file. */
void ReadMeshFromGmshFile( std::string gmshFileName )
{
	std::stringstream msg("");
	
	std::ifstream gmshFileInput(gmshFileName.c_str());
	//open file for reading
	if (!gmshFileInput){
		msg.str(" ");
		msg << "Cannot open Gmsh file for reading."<<std::endl <<
			"Please check the filename and try again."<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::exit(1);
	}
	char str[255];
	gmshFileInput.getline(str,255);
	
	// read in the nodes /
	while ( strcmp(str,"$Nodes") ){
		gmshFileInput.getline(str,255); // skip all lines until we get to "$Nodes"
	}
	unsigned int numberOfNodes;
	gmshFileInput >> numberOfNodes; // the next line will be the number of nodes
	
	// create the data image
	DataImagePointsPointer		points				= DataImagePointsPointer::New();
	DataImagePointer			meshImage			= DataImagePointer::New();
	DataImagePixelPointer		displacementData	= DataImagePixelPointer::New();
	DataImagePixelPointer		optimizerData		= DataImagePixelPointer::New();
	
	meshImage->SetPoints( points );
	
	displacementData->SetNumberOfComponents(3);
	displacementData->SetName("Displacement");
	meshImage->GetPointData()->AddArray( displacementData );
	
	optimizerData->SetNumberOfComponents(1);
	optimizerData->SetName("Optimizer Value");
	meshImage->GetPointData()->AddArray( optimizerData );
	
	meshImage->GetPoints()->SetNumberOfPoints( (vtkIdType)numberOfNodes );
	for (unsigned int i = 0; i < numberOfNodes; ++i){ // run through the number of nodes and put them in the pointset
		unsigned int pointNumber;
		float x, y, z;
		gmshFileInput >> pointNumber >> x >> y >> z;
		meshImage->GetPoints()->SetPoint(i,x,y,z);
	}
	
	// read in the elements. /
	while ( strcmp(str,"$Elements") ){
		gmshFileInput.getline(str,255); // skip all lines until we get to "$Elements"
	}
	unsigned int numberOfElements;
	gmshFileInput >> numberOfElements; // the next line will be the number of elements
	
	for ( unsigned int i = 0; i < numberOfElements; ++i ){
		unsigned int elNo, elTyp;
		gmshFileInput >> elNo >> elTyp;
		
		//~ if (elTyp == 1){ // 2 node line
			//~ unsigned int elNTags;
			//~ gmshFileInput >> elNTags;
			//~ for (unsigned int j = 0; j < elNTags+1; ++j){
				//~ gmshFileInput.ignore(255,' ');
			//~ }
			//~ int nodes[2];
			//~ gmshFileInput >> nodes[0] >> nodes[1];
			//~ vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1 };
			//~ meshImage->InsertNextCell(3, 2, ptIds);
		//~ }
		
		//~ if (elTyp == 2){ // 3 node triangle
			//~ unsigned int elNTags;
			//~ gmshFileInput >> elNTags;
			//~ for (unsigned int j = 0; j < elNTags+1; ++j){
				//~ gmshFileInput.ignore(255,' ');
			//~ }
			//~ int nodes[3];
			//~ gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2];
			//~ vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1 };
			//~ meshImage->InsertNextCell(5, 3, ptIds);
		//~ }
		
		if (elTyp == 4){ // linear tet element
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0;j < elNTags+1; ++j){ // ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[4];
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3]; 
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1 };
			meshImage->InsertNextCell(10, 4, ptIds);
		}
		
		//~ if (elTyp == 8){ // quadratic line element
			//~ std::stringstream msg("");
			//~ msg << "Warning: Qadratic line element in gmsh file."<<std::endl << "Qadratic lines are not supported" <<std::endl <<
				//~ "Continuing with the file input."<<std::endl; // not supported because they give bad jaboians in the strain calcs
			//~ this->WriteToLogfile( msg.str() );
			//~ gmshFileInput.ignore(255,'\n');
			//~ continue;
		//~ }
		
		//~ if (elTyp == 9){ //quadratic triangle element
			//~ unsigned int elNTags;
			//~ gmshFileInput >> elNTags;
			//~ for (unsigned int j = 0; j< elNTags+1; ++j){ //ignore tags
				//~ gmshFileInput.ignore(255,' ');
			//~ }
			//~ int nodes[6];
			//~ gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5];
			//~ vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1, nodes[4]-1, nodes[5]-1};
			//~ meshImage->InsertNextCell(22,6,ptIds);
		//~ }
		
		if (elTyp == 11){ // quatratic tet element
			unsigned int elNTags;
			gmshFileInput >> elNTags;
			for (unsigned int j = 0;j < elNTags+1; ++j){ // ignore tags
				gmshFileInput.ignore(255,' ');
			}
			int nodes[10];//, node1, node2, node3, node4, node5, node6, node7, node8, node9;
			gmshFileInput >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7] >> nodes[9] >> nodes[8]; // nodes[8] and modes[9] intensionally switched due to vtk/gmsh node numbering differences
			vtkIdType ptIds[] = { nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1, nodes[4]-1, nodes[5]-1, nodes[6]-1, nodes[7]-1, nodes[8]-1, nodes[9]-1 };
			meshImage->InsertNextCell(24, 10, ptIds);
		}
		
		//~ if (elTyp == 15){ // 3D vertex
			//~ unsigned int elNTags;
			//~ gmshFileInput >> elNTags;
			//~ for (unsigned int j = 0; j< elNTags+1; ++j){ //ignore tags
				//~ gmshFileInput.ignore(255,' ');
			//~ }
			//~ int nodes[1];
			//~ gmshFileInput >> nodes[0];
			//~ vtkIdType ptIds[] = { nodes[0]-1 };
			//~ meshImage->InsertNextCell(1,1,ptIds);
		//~ }
		
		//~ if (elTyp != 1 && elTyp != 2 && elTyp != 4 && elTyp != 8 && elTyp != 9 && elTyp != 11 && elTyp != 15){ // only support the above element types
			//~ std::stringstream msg("");
			//~ msg << "Unsupported element types detected!"<<std::endl << "Unsupported type = "<<elTyp<<std::endl <<
				//~ "Continuing with the file input."<<std::endl;
			//~ this->WriteToLogfile( msg.str() );
			//~ gmshFileInput.ignore(255,'\n');
			//~ continue;
		//~ }
		
		if (elTyp != 4 && elTyp != 11){ // only support the above element types
			std::stringstream msg("");
			msg << "Unsupported element types detected!"<<std::endl << "Unsupported type = "<<elTyp<<std::endl <<
				"Continuing with the file input."<<std::endl;
			this->WriteToLogfile( msg.str() );
			gmshFileInput.ignore(255,'\n');
			continue;
		}	
	}
		
	for ( unsigned int i = 0; i < numberOfNodes; ++i ){
		meshImage->GetPointData()->GetArray("Displacement")->InsertNextTuple3( 0, 0, 0 );
		meshImage->GetPointData()->GetArray("Optimizer Value")->InsertNextTuple1( 0 );
	}
	
	this->SetDataImage( meshImage );
}

/** A function to fill a mesh with an single value. */
void SetMeshInitialValue( double initialData[] )
{
	for ( unsigned int i = 0; i < this->m_DataImage->GetNumberOfPoints(); ++i ){
		this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i , initialData );
	}

}

/** Get a pixel value from the mesh using an index */
void GetMeshPixelValueFromIndex( vtkIdType index, double *pixel )
{
	this->m_DataImage->GetPointData()->GetArray("Displacement")->GetTuple( index, pixel );
}

/** Get the value of the optimizer at a certain point. */
void GetMeshPixelOptimizerFromIndex( vtkIdType index, double *opt )
{
	this->m_DataImage->GetPointData()->GetArray("Optimizer Value")->GetTuple( index, opt );
}

/** Set a pixel value in the mesh using an index */
void SetMeshPixelValueFromIndex(vtkIdType index, double *pixel )
{
	this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( index, pixel );
}

/** Set the value of the optimizer at a certain point. */
void SetMeshPixelOptimizerFromIndex( vtkIdType index, double *opt )
{
	this->m_DataImage->GetPointData()->GetArray("Optimizer Value")->SetTuple( index, opt );
}

/** Get a point by index from the mesh. */
void GetMeshPointLocationFromIndex( vtkIdType index, double  *point )
{
	this->m_DataImage->GetPoint(index,point);
}

/** Set the data Image. */
void SetDataImage( vtkUnstructuredGrid *initialDataImage )
{
	if (this->m_DataImage.GetPointer() != initialDataImage){
		this->m_DataImage = initialDataImage;
	}
	this->KDTreeSetAndBuild();
}

/** Get the pointer to the data image. */
DataImagePointer GetDataImage()
{
	return this->m_DataImage;
}

/** a function to execute the analysis. */
void ExecuteDIC()
{
	std::stringstream msg(""); // test if everything is set
	if( !this->m_Registration ){
		msg.str("");
		msg << "Registration not set.  Please define and set the\nregistration using the SetRegistrationMethod() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	if( !this->m_FixedImage ){
		msg.str("");
		msg << "Fixed image not set.  Please define and set the\n fixed image using the SetFixedImage() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	if( !this->m_MovingImage ){
		msg.str("");
		msg << "Moving image not set.  Please define and set the\n moving image using the SetMovingImage() method.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();		
	}
	
	if( !this->m_DataImage ){
		msg.str("");
		msg << "The data image not set.  Please define and set the\n initial data image using either the ReadMeshFromGmshFile or\n SetDataImage() methods.\n";
		this->WriteToLogfile( msg.str() );
		std::abort();		
	}
	
	if( this->m_FixedImageRegionList.empty() ){ // if the region list is empty, create full region lists
		this->CalculateInitialFixedImageRegionList();
	}
	if( this->m_MovingImageRegionList.empty() ){
		this->CalculateInitialMovingImageRegionList();
	}
		
	struct tm * timeValue; // record the time
	std::time_t rawTime;
	std::time( &rawTime );
	timeValue = std::localtime( &rawTime );
	msg.str("");
	msg << "Starting DIC at: "<<std::asctime( timeValue )<<std::endl;
	this->WriteToLogfile( msg.str() );
	
	// visit every point in the points list
	unsigned int nMeshPoints = this->m_pointsList->GetNumberOfIds();
	for( unsigned int i = 0; i<nMeshPoints; ++i){
		vtkIdType pointId = this->m_pointsList->GetId( i );
		
		msg.str("");
		msg << "Starting image registraion for point: "<<i+1<<" of "<<nMeshPoints<<" (mesh index "<<pointId<<")"<<std::endl;
		this->WriteToLogfile( msg.str() );
		
		std::time_t rawTime; // record the time for each DVC
		struct tm * timeValue;
		std::time( &rawTime );
		timeValue = std::localtime( &rawTime );
		msg.str("");
		msg << "Time: "<<std::asctime( timeValue );
		this->WriteToLogfile( msg.str() );
		
		FixedImageRegionType	*fixedRegion = new FixedImageRegionType ; // get the fixed image from the fixed image list
		FixedImagePointer		fixedImage;
		fixedRegion = this->GetFixedImageRegionFromIndex( i );
		fixedImage = this->GetFixedROIAsImage( fixedRegion );
		this->SetFixedROIImage( fixedImage );
		
		MovingImageRegionType	*movingRegion = new MovingImageRegionType; // get the moving image from the moving image list
		MovingImagePointer		movingImage;
		movingRegion = this->GetMovingImageRegionFromIndex( i );
		movingImage = this->GetMovingROIAsImage( movingRegion );
		this->SetMovingROIImage( movingImage );
		
		this->SetTransformToIdentity(); // Set the transform to do nothing
		
		// initialize the transform to perform rotations about the fixed region center
		this->m_TransformInitializer->SetFixedImage( this->m_Registration->GetFixedImage() ); 
		this->m_TransformInitializer->SetMovingImage( this->m_Registration->GetMovingImage() );
		this->m_TransformInitializer->SetTransform( this->m_Transform );
		this->m_TransformInitializer->GeometryOn();
		this->m_TransformInitializer->InitializeTransform();
		this->m_Registration->SetInitialTransformParameters( this->m_Transform->GetParameters() );
		
		double	*displacementData = new double[3];  // set the initial displacement
		this->GetMeshPixelValueFromIndex( pointId, displacementData );
		this->SetInitialDisplacement( displacementData );
		
		msg.str("");
		msg <<"Current transform: "<<this->m_Registration->GetInitialTransformParameters()<<std::endl;
		this->WriteToLogfile( msg.str() );
		
		// if the optimizer is the lbfgsb then set the bounds based on teh current displacement
		if( !strcmp(this->m_Registration->GetOptimizer()->GetNameOfClass(),"LBFGSBOptimizer") ){
			unsigned int nParameters = this->m_Registration->GetTransform()->GetNumberOfParameters();
			itk::Array< long > boundSelect( nParameters );
			itk::Array< double > upperBound( nParameters );
			itk::Array< double > lowerBound( nParameters );
			
			boundSelect.Fill( 0 );
			boundSelect[nParameters-3] = 2;
			boundSelect[nParameters-2] = 2;
			boundSelect[nParameters-1] = 2;
			typename FixedImageType::SpacingType imageSpacing = this->m_Registration->GetFixedImage()->GetSpacing();
			upperBound[nParameters-3] = *displacementData + .5*imageSpacing[0];
			upperBound[nParameters-2] = *(displacementData+1) + .5*imageSpacing[1];
			upperBound[nParameters-1] = *(displacementData+2) + .5*imageSpacing[2];
			lowerBound.Fill( 0 );
			lowerBound[nParameters-3] = *displacementData - .5*imageSpacing[0];
			lowerBound[nParameters-2] = *(displacementData+1) - .5*imageSpacing[1];
			lowerBound[nParameters-1] = *(displacementData+2) - .5*imageSpacing[2];
			
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetBoundSelection( boundSelect );
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetUpperBound( upperBound );
			reinterpret_cast<itk::LBFGSBOptimizer *>(this->m_Registration->GetOptimizer())->SetLowerBound( lowerBound );
		}
		
		// update the registration
		this->UpdateRegionRegistration();
		
		// output the results
		double *lastDisp = new double[3];
		this->GetLastDisplacement( lastDisp );
		this->SetMeshPixelValueFromIndex( pointId, lastDisp );
		double *lastOpt = new double;
		this->GetLastOptimizer( lastOpt );
		this->SetMeshPixelOptimizerFromIndex( pointId, lastOpt );
		msg.str("");
		msg << "Final displacement value: ("<<*lastDisp<<", "<<*(lastDisp+1)<<", "<<*(lastDisp+2)<<")"<<std::endl <<
			"Optimizer stop condition: " << this->m_Registration->GetOptimizer()->GetStopConditionDescription() << std::endl <<
			"Final optimizer value: "<<*lastOpt<<std::endl<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::string debugFile = this->m_OutputDirectory + "/debug.vtk";
		this->WriteMeshToVTKFile( debugFile );
	}
}

/** A function to write the mesh data to a VTK ASCII file. */
void WriteMeshToVTKFile(std::string outFile)
{
	vtkSmartPointer<vtkUnstructuredGridWriter>	writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName( outFile.c_str() );
	writer->SetFileTypeToASCII();
	writer->SetInput( this->m_DataImage );
	writer->Update();
}

/** A function to build the KDTree locator from the data image.  This method
 * is called when the data image is set. */
void KDTreeSetAndBuild()
{
	this->m_KDTree->SetDataSet( this->m_DataImage );
	this->m_KDTree->BuildLocatorFromPoints( this->m_DataImage->GetPoints() );
}

/** A function to find the values that are outside a given bounds 
 * compared to their connected neighbours. */
void CreateNewRegionListFromBadPixels()
{	
	this->m_FixedImageRegionList.clear(); // clear the regions for calculation
	this->m_MovingImageRegionList.clear();
	this->m_pointsList->Reset(); // clear the point list
	
	DataImagePointer tempImage = DataImagePointer::New(); // make a copy of the image to work with
	tempImage->DeepCopy(this->m_DataImage);
	
	// visit every point in the image
	unsigned int numberOfPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < numberOfPoints; ++i ){
		
		// find connected points
		vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
		this->GetNeighbouringPoints( i, pointList );

		// calc the average and standard deviabtions of the displacements stored in the neighbourhood points
		double averageValue[3] = {0,0,0};
		double averageMag = 0;
		double stDevValue[3] = {0,0,0};
		double stDevMag = 0;
		this->GetDisplacementStats(i, pointList, averageValue, averageMag, stDevValue, stDevMag, tempImage );
		
		// get the displacement value stored in the current point
		double *currentValue = new double[3];
		currentValue = tempImage->GetPointData()->GetArray("Displacement")->GetTuple( i );
	
		// check if there the pixel is valid - if any component is out of value it is considered invalid
		if( !DisplacementValid(currentValue, averageValue, averageMag, stDevValue, stDevMag ) ){
			
			double *averagedDisplacement = new double[3];
			averagedDisplacement = CalculateDisplacementWeightedMovingAverage( 2, 0, i, tempImage );
			this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i, averagedDisplacement );  // apply a smoothed value to the mesh to be used in the next evaluation		
				
			double *movingImageCenterLocation = new double[3];
			movingImageCenterLocation = this->GetMovingImageRegionLocationFromIndex( i ); // new moving image centre
			
			MovingImageRegionType *currentMovingRegion = new MovingImageRegionType;
			this->GetMovingImageRegionFromLocation( currentMovingRegion, movingImageCenterLocation );
			this->PushRegionOntoMovingImageRegionList( currentMovingRegion ); // new moving image
			
			FixedImageRegionType *currentFixedRegion = new FixedImageRegionType;  // calculated the new fixed region
			this->GetFixedImageRegionFromLocation( currentFixedRegion, this->m_DataImage->GetPoint( i ) );
			this->PushRegionOntoFixedImageRegionList( currentFixedRegion );
			
			this->m_pointsList->InsertNextId( i ); // add the current point to the next round of evaluation
		}
	}
}

/** A function to calcluate the component and mangnitued averages of point
 * values from a list of point IDs. 
 * rPoint = reference point
 * points = list of points for stats calculation
 * vectorAverage = average value of each component (eg, vectorAverage[0] = average(x_displacement) 
 * magAverage = average magnitude of vectorAverage
 * vectorStDev = standard deviation of each component
 * magStDev = magnitude of the standard deviations
 * image = pointer to an unstructured grid image */ 
void GetDisplacementStats(vtkIdType rPoint, vtkSmartPointer<vtkIdList> points, double vectorAverage[3], double &magAverage, double vectorStDev[3], double &magStDev, DataImagePointer image)
{
	// visit each point in the list and calculate the average
	unsigned int nPoints = points->GetNumberOfIds();
	int denominator = 0;
	for( unsigned int i = 0; i < nPoints; ++i){
		
		vtkIdType pointId = points->GetId( i );
		if ( pointId == rPoint) continue; // don't include the current point in the average.
		denominator++; // increment the denominator only once passed the current point check
		
		double *cPoint = new double[3];
		image->GetPointData()->GetArray("Displacement")->GetTuple( pointId, cPoint ); // get the current point value
		
		vectorAverage[0] = vectorAverage[0] + cPoint[0];
		vectorAverage[1] = vectorAverage[1] + cPoint[1];
		vectorAverage[2] = vectorAverage[2] + cPoint[2];
		
		magAverage = magAverage + std::sqrt( std::pow(cPoint[0],2) + std::pow(cPoint[1],2) + std::pow(cPoint[2],2) );
	}
	
	vectorAverage[0] = vectorAverage[0]/denominator;
	vectorAverage[1] = vectorAverage[1]/denominator;
	vectorAverage[2] = vectorAverage[2]/denominator; 
	magAverage = magAverage/denominator;
	
	double vectDiff[3] = {0, 0, 0};
	double magDiff = 0;
	denominator = 0;
	for( unsigned int i = 0; i <nPoints; ++i){
		
		vtkIdType pointId = points->GetId( i );
		if ( pointId == rPoint) continue; // don't include the current point in the average.
		denominator++; // increment the denominator only once passed the current point check
		
		double *cPoint = new double[3];
		image->GetPointData()->GetArray("Displacement")->GetTuple( pointId, cPoint ); // get the current point value
		
		double cMag = std::sqrt( std::pow(cPoint[0],2) + std::pow(cPoint[1],2) + std::pow(cPoint[2],2) );
		
		// the sum of all the squared distances from the mean
		vectDiff[0] = vectDiff[0] + std::pow( (cPoint[0]-vectorAverage[0]), 2);
		vectDiff[1] = vectDiff[1] + std::pow( (cPoint[1]-vectorAverage[1]), 2);
		vectDiff[2] = vectDiff[2] + std::pow( (cPoint[2]-vectorAverage[2]), 2);
		
		magDiff = magDiff + std::pow( (cMag - magAverage), 2);
	}
	
	vectorStDev[0] = std::sqrt(vectDiff[0]/denominator);
	vectorStDev[1] = std::sqrt(vectDiff[1]/denominator);
	vectorStDev[2] = std::sqrt(vectDiff[2]/denominator);
	
	magStDev = std::sqrt(magDiff/denominator);
}

void ReplaceBadDisplacementPixels( double sigma, double mean, vtkIdList *replacedList)
{
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	
	// empty the list of replaced points
	replacedList->Reset();
	
	// visit every point in the image
	unsigned int numberOfPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < numberOfPoints; ++i ){

		// get connected points
		vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();		
		this->GetNeighbouringPoints( i, pointList );

		// calc the average and stDev of the points in the list
		double vectorAverage[3] = { 0 };
		double magAverage = 0;
		double vectorStDev[3] = { 0 };
		double magStDev = 0;		
		this->GetDisplacementStats( i, pointList, vectorAverage, magAverage, vectorStDev, magStDev, tempImage );
			
		// get the value of the current point
		double *currentValue = new double[3];
		tempImage->GetPointData()->GetArray("Displacement")->GetTuple( i, currentValue );
		
		// check if there the pixel is valid - if any component is out of value it is considered invalid
		if( !this->DisplacementValid( currentValue, vectorAverage, magAverage, vectorStDev, magStDev ) ){	
			double *newDisplacement = new double[3];
			newDisplacement = CalculateDisplacementWeightedMovingAverage( sigma, mean, i, tempImage );
			this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i, newDisplacement );  // apply pixel value to the mesh
			replacedList->InsertNextId( i );
		}
	}
}

bool DisplacementValid(double value[3], double vectorAverage[3], double &magAverage, double vectorStDev[3], double &magStDev)
{
		if (
		std::fabs(value[0] - vectorAverage[0]) > this->m_displacementErrorTolerance * vectorStDev[0] ||
		std::fabs(value[1] - vectorAverage[1]) > this->m_displacementErrorTolerance * vectorStDev[1] ||
		std::fabs(value[2] - vectorAverage[2]) > this->m_displacementErrorTolerance * vectorStDev[2] )
		{
			return false;
		}
		else
		{
			return true;
		}
	
}

void ReplaceBadStrainPixels( double sigma, double mean, vtkIdList *replacedList)
{
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	
	// empty the list of replaced points
	replacedList->Reset(); 
	
	// visit every point in the image
	unsigned int numberOfPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < numberOfPoints; ++i ){

		// get connected points
		vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();		
		this->GetNeighbouringPoints( i, pointList );

		// calc the average and stDev of the points in the list
		double averageValue[9] = { 0 };
		double stDevValue[9] = { 0 };
		this->GetStrainStats(i, pointList, averageValue, stDevValue );
		
		// get the value of the current point
		double *currentValue = new double[9];
		currentValue = tempImage->GetPointData()->GetArray("Strain")->GetTuple( i );
		
		// check if there the pixel is valid - if any component is out of value it is considered invalid
		if( !this->StrainValid(currentValue, averageValue, stDevValue) ){	
			double *newStrain = new double[9];
			newStrain = CalculateStrainWeightedMovingAverage( sigma, mean, i, tempImage );
			this->m_DataImage->GetPointData()->GetArray("Strain")->SetTuple( i, newStrain );  // apply pixel value to the mesh
			replacedList->InsertNextId( i );
		}
	}
}

void GetStrainStats( vtkIdType c_point, vtkSmartPointer<vtkIdList> points, double vectorAverage[9], double vectorStDev[9] )
{
	unsigned int nPoints = points->GetNumberOfIds();
	int denominator = 0;
	for( unsigned int i = 0; i < nPoints; ++i){
		vtkIdType pointId = points->GetId( i );
		if ( pointId == c_point) continue; // don't include the current point in the average.
		denominator++;
		double *cPoint = new double[9];
		this->m_DataImage->GetPointData()->GetArray("Strain")->GetTuple( pointId, cPoint );
		
		vectorAverage[0] = vectorAverage[0] + cPoint[0];
		vectorAverage[1] = vectorAverage[1] + cPoint[1];
		vectorAverage[2] = vectorAverage[2] + cPoint[2];
		vectorAverage[3] = vectorAverage[3] + cPoint[3];
		vectorAverage[4] = vectorAverage[4] + cPoint[4];
		vectorAverage[5] = vectorAverage[5] + cPoint[5];
		vectorAverage[6] = vectorAverage[6] + cPoint[6];
		vectorAverage[7] = vectorAverage[7] + cPoint[7];
		vectorAverage[8] = vectorAverage[8] + cPoint[8];
	}
	
	vectorAverage[0] = vectorAverage[0]/denominator;
	vectorAverage[1] = vectorAverage[1]/denominator;
	vectorAverage[2] = vectorAverage[2]/denominator;
	vectorAverage[3] = vectorAverage[3]/denominator;
	vectorAverage[4] = vectorAverage[4]/denominator;
	vectorAverage[5] = vectorAverage[5]/denominator;
	vectorAverage[6] = vectorAverage[6]/denominator;
	vectorAverage[7] = vectorAverage[7]/denominator;
	vectorAverage[8] = vectorAverage[8]/denominator;
	
	double vectDiff[9] = { 0 };
	denominator = 0;
	for( unsigned int i = 0; i <nPoints; ++i){
		vtkIdType pointId = points->GetId( i );
		if ( pointId == c_point) continue; // don't include the current point in the average.
		denominator++;
		double *cPoint = new double[9];
		this->m_DataImage->GetPointData()->GetArray("Strain")->GetTuple( pointId, cPoint );
		
		vectDiff[0] = vectDiff[0] + std::pow( (cPoint[0] - vectorAverage[0]), 2);
		vectDiff[1] = vectDiff[1] + std::pow( (cPoint[1] - vectorAverage[1]), 2);
		vectDiff[2] = vectDiff[2] + std::pow( (cPoint[2] - vectorAverage[2]), 2);
		vectDiff[3] = vectDiff[3] + std::pow( (cPoint[3] - vectorAverage[3]), 2);
		vectDiff[4] = vectDiff[4] + std::pow( (cPoint[4] - vectorAverage[4]), 2);
		vectDiff[5] = vectDiff[5] + std::pow( (cPoint[5] - vectorAverage[5]), 2);
		vectDiff[6] = vectDiff[6] + std::pow( (cPoint[6] - vectorAverage[6]), 2);
		vectDiff[7] = vectDiff[7] + std::pow( (cPoint[7] - vectorAverage[7]), 2);
		vectDiff[8] = vectDiff[8] + std::pow( (cPoint[8] - vectorAverage[8]), 2);
	}
	
	vectorStDev[0] = std::sqrt( vectDiff[0]/denominator );
	vectorStDev[1] = std::sqrt( vectDiff[1]/denominator );
	vectorStDev[2] = std::sqrt( vectDiff[2]/denominator );
	vectorStDev[3] = std::sqrt( vectDiff[3]/denominator );
	vectorStDev[4] = std::sqrt( vectDiff[4]/denominator );
	vectorStDev[5] = std::sqrt( vectDiff[5]/denominator );
	vectorStDev[6] = std::sqrt( vectDiff[6]/denominator );
	vectorStDev[7] = std::sqrt( vectDiff[7]/denominator );
	vectorStDev[8] = std::sqrt( vectDiff[8]/denominator );
}

bool StrainValid(double currentValue[9], double averageValue[9], double stDevValue[9])
{
	if (std::fabs(currentValue[0] - averageValue[0]) > this->m_strainErrorTolerance * stDevValue[0] ||
		std::fabs(currentValue[1] - averageValue[1]) > this->m_strainErrorTolerance * stDevValue[1] ||
		std::fabs(currentValue[2] - averageValue[2]) > this->m_strainErrorTolerance * stDevValue[2] ||
		std::fabs(currentValue[3] - averageValue[3]) > this->m_strainErrorTolerance * stDevValue[3] ||
		std::fabs(currentValue[4] - averageValue[4]) > this->m_strainErrorTolerance * stDevValue[4] ||
		std::fabs(currentValue[5] - averageValue[5]) > this->m_strainErrorTolerance * stDevValue[5] ||
		std::fabs(currentValue[6] - averageValue[6]) > this->m_strainErrorTolerance * stDevValue[6] ||
		std::fabs(currentValue[7] - averageValue[7]) > this->m_strainErrorTolerance * stDevValue[7] ||
		std::fabs(currentValue[8] - averageValue[8]) > this->m_strainErrorTolerance * stDevValue[8] )
		{
			return false;
		}
		else
		{
			return true;
		}
}

/** A function to calculate the strains. */
void GetStrains()
{
	// calculate derivatives
	vtkCellDerivatives *derivativeCalculator = vtkCellDerivatives::New();
	derivativeCalculator->SetInput( this->m_DataImage );
	derivativeCalculator->SetTensorModeToComputeStrain();
	derivativeCalculator->SetVectorModeToComputeGradient();
	
	derivativeCalculator->SetInputArrayToProcess(1,0,0,0,"Displacement"); // the CellDerivatives calculator sets the array idx to 0 for scalars and 1 for vector inputs. We want it to operate on vectors, so idx is set to 1.  Further, giving it the name of the array allows it to verify that it is operating on the correct array.
	derivativeCalculator->Update();
	this->SetDataImage( derivativeCalculator->GetUnstructuredGridOutput() );
	
	// move the data from cells to nodes
	vtkCellDataToPointData	*dataTransformer = vtkCellDataToPointData::New();
	dataTransformer->SetInput( this->m_DataImage );
	dataTransformer->PassCellDataOn();
	dataTransformer->SetInputArrayToProcess(0,0,0,1,"Strain");// Strain is on array zero of the cells in m_DataImage. The 4th is an enum to tell the algorithm the data is stored in cells, the last in the name of the array.
	dataTransformer->Update();
	
	this->SetDataImage( dataTransformer->GetUnstructuredGridOutput() );
}

/** A function to calculate the principal strains. */
void GetPrincipalStrains()
{
	vtkSmartPointer<vtkMath>		mathAlgorithm = vtkSmartPointer<vtkMath>::New();
	// First find the principal strains for the point data
	vtkIdType nPoints = this->m_DataImage->GetPointData()->GetArray("Strain")->GetNumberOfTuples();
	
	DataImagePixelPointer V0;
	DataImagePixelPointer V1;
	DataImagePixelPointer V2;
	DataImagePixelPointer val0;
	DataImagePixelPointer val1;
	DataImagePixelPointer val2;
	
	// create the eigenvector containers if they don't exist
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 1") ){
		V0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 1" ) );

	}
	else{
		V0 = DataImagePixelPointer::New();
		V0->SetNumberOfComponents(3);
		V0->SetName("Principal Strain Vector 1");
		V0->Allocate(nPoints);
		V0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V0 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 2") ){
		V1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 2" ) );
	}
	else{
		V1 = DataImagePixelPointer::New();
		V1->SetNumberOfComponents(3);
		V1->SetName("Principal Strain Vector 2");
		V1->Allocate(nPoints);
		V1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V1 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Vector 3") ){
		V2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Vector 3" ) );
	}
	else{
		V2 = DataImagePixelPointer::New();
		V2->SetNumberOfComponents(3);
		V2->SetName("Principal Strain Vector 3");
		V2->Allocate(nPoints);
		V2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( V2 );
	}

	// create the eigenvalue containers if they don't exist
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 1") ){
		val0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray( "Principal Strain Value 1" ) );
	}
	else{
		val0 = DataImagePixelPointer::New();
		val0->SetNumberOfComponents(1);
		val0->SetName("Principal Strain Value 1");
		val0->Allocate(nPoints);
		val0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val0 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 2") ){
		val1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 2") );
	}
	else{
		val1 = DataImagePixelPointer::New();
		val1->SetNumberOfComponents(1);
		val1->SetName("Principal Strain Value 2");
		val1->Allocate(nPoints);
		val1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val1 );
	}
	if ( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 3") ){
		val2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetPointData()->GetArray("Principal Strain Value 3") );
	}
	else{
		val2 = DataImagePixelPointer::New();
		val2->SetNumberOfComponents(1);
		val2->SetName("Principal Strain Value 3");
		val2->Allocate(nPoints);
		val2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetPointData()->AddArray( val2 );
	}	
	
	for ( unsigned int i = 0; i < nPoints; ++i ){
		// Get the point tensor
		double *cTensorRaw = this->m_DataImage->GetPointData()->GetArray("Strain")->GetTuple( i );
		// reform tensor into an appropreate array
		double *cTensor[3];
		double cT0[3];
		double cT1[3];
		double cT2[3];

		cTensor[0] = cT0; cTensor[1] = cT1; cTensor[2] = cT2;

		cT0[0] = *cTensorRaw;
		cT0[1] = *(cTensorRaw+1);
		cT0[2] = *(cTensorRaw+2);	
		cT1[0] = *(cTensorRaw+3);
		cT1[1] = *(cTensorRaw+4);
		cT1[2] = *(cTensorRaw+5);
		cT2[0] = *(cTensorRaw+6);
		cT2[1] = *(cTensorRaw+7);
		cT2[2] = *(cTensorRaw+8);
				
		// Get the eigen-values and -vectors
		double eigenValues[3]; // The Eigenvalues
		double *eigenVectors[3]; // vector of pointers to eigenvectors
		double v0[3]; // eigen vector 1
		double v1[3]; // eigen vector 2
		double v2[3]; // eigen vector 3
		eigenVectors[0] = v0; eigenVectors[1] = v1; eigenVectors[2] = v2;
		
		// perform the calculation of the principal strains
		mathAlgorithm->Jacobi(cTensor,eigenValues,eigenVectors);
		
		// save the values and vectors into the data image
		V0->SetTuple( i, eigenVectors[0] );
		val0->SetTuple( i, &eigenValues[0] );
		V1->SetTuple( i, eigenVectors[1] );
		val1->SetTuple( i, &eigenValues[1] );
		V2->SetTuple( i, eigenVectors[2] );
		val2->SetTuple( i, &eigenValues[2] );
	}
	
	// next find the principal strains for the cell data
	nPoints = this->m_DataImage->GetCellData()->GetArray("Strain")->GetNumberOfTuples();
	
	// create the eigenvector containers if they don't exist
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 1") ){
		V0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 1") );
	}
	else{
		V0 = DataImagePixelPointer::New();
		V0->SetNumberOfComponents(3);
		V0->SetName("Principal Strain Vector 1");
		V0->Allocate(nPoints);
		V0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V0 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 2") ){
		V1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 2") );
	}
	else{
		V1 = DataImagePixelPointer::New();
		V1->SetNumberOfComponents(3);
		V1->SetName("Principal Strain Vector 2");
		V1->Allocate(nPoints);
		V1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V1 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 3") ){
		V2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Vector 3") );
	}
	else{
		V2 = DataImagePixelPointer::New();
		V2->SetNumberOfComponents(3);
		V2->SetName("Principal Strain Vector 3");
		V2->Allocate(nPoints);
		V2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( V2 );
	}

	// create the eigenvalue containers if they don't exist
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 1") ){
		val0 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 1") );
	}
	else{
		val0 = DataImagePixelPointer::New();
		val0->SetNumberOfComponents(1);
		val0->SetName("Principal Strain Value 1");
		val0->Allocate(nPoints);
		val0->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val0 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 2") ){
		val1 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 2") );
	}
	else{
		val1 = DataImagePixelPointer::New();
		val1->SetNumberOfComponents(1);
		val1->SetName("Principal Strain Value 2");
		val1->Allocate(nPoints);
		val1->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val1 );
	}
	if ( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 3") ){
		val2 = vtkDoubleArray::SafeDownCast( this->m_DataImage->GetCellData()->GetArray("Principal Strain Value 3") );
	}
	else{
		val2 = DataImagePixelPointer::New();
		val2->SetNumberOfComponents(1);
		val2->SetName("Principal Strain Value 3");
		val2->Allocate(nPoints);
		val2->SetNumberOfTuples(nPoints);
		this->m_DataImage->GetCellData()->AddArray( val2 );
	}
	
	for ( unsigned int i = 0; i < nPoints; ++i ){
		// Get the point tensor
		double *cTensorRaw = this->m_DataImage->GetCellData()->GetArray("Strain")->GetTuple( i );

		// reform tensor into an appropreate array
		double *cTensor[3];
		double cT0[3];
		double cT1[3];
		double cT2[3];

		cTensor[0] = cT0; cTensor[1] = cT1; cTensor[2] = cT2;

		cT0[0] = *cTensorRaw;
		cT0[1] = *(cTensorRaw+1);
		cT0[2] = *(cTensorRaw+2);	
		cT1[0] = *(cTensorRaw+3);
		cT1[1] = *(cTensorRaw+4);
		cT1[2] = *(cTensorRaw+5);
		cT2[0] = *(cTensorRaw+6);
		cT2[1] = *(cTensorRaw+7);
		cT2[2] = *(cTensorRaw+8);

		// Get the eigen-values and -vectors
		double eigenValues[3]; // The Eigenvalues
		double *eigenVectors[3]; // vector of pointers to eigenvectors
		double v0[3]; // eigen vector 1
		double v1[3]; // eigen vector 2
		double v2[3]; // eigen vector 3
		eigenVectors[0] = v0; eigenVectors[1] = v1; eigenVectors[2] = v2;

		// perform the calculation of the principal strains
		mathAlgorithm->Jacobi(cTensor,eigenValues,eigenVectors);

		// save the values and vectors into the data image
		V0->SetTuple( i, eigenVectors[0] );
		val0->SetTuple( i, &eigenValues[0] );
		V1->SetTuple( i, eigenVectors[1] );
		val1->SetTuple( i, &eigenValues[1] );
		V2->SetTuple( i, eigenVectors[2] );
		val2->SetTuple( i, &eigenValues[2] );
	}
	
}

/** A function to smooth the mesh image.  This method takes the N closest 
 * points (includeing the point in question) and takes the point with the
 * median magnitude as the value. N must be odd or 1 will be added to it.
 * NOTE: This filter has been replaced by "ReplaceBadDisplacementPixels"
 * and "ReplaceBadStrainPixels." It has not been updated or debugged - 
 * USE AT YOUR OWN RISK (or update and debug). */
void MedianImageFilter(unsigned int N)
{
	// if N is even, add one
	N = (N % 2 == 1) ? N += 1 : N;
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	// loop through the points of both the images
	unsigned int nPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < nPoints; ++i ){
		double *cPoint = tempImage->GetPoint( i );
		// find the 7 closest points to the position of each point
		vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New();
		//~ this->m_KDTree->FindClosestNPoints( N, cPoint, pointList );
		this->GetNeighbouringPoints( i, pointList );
		// calculate magnitudes
		double neighbourhoodMagnitudes[N];
		std::vector<double> sortingArray; // array for magnitude sorting
		for( unsigned int j = 0; j < N; ++j){
			vtkIdType cId = pointList->GetId( j );
			double *tempPixel = tempImage->GetPointData()->GetArray("Displacement")->GetTuple( cId );
			double cMag = sqrt( pow(*tempPixel,2) + pow(*(tempPixel+1),2) + pow(*(tempPixel+2),2) );
			neighbourhoodMagnitudes[j] = cMag ;
			sortingArray.push_back( cMag );
		}
		// sort the magnitudes
		std::sort( sortingArray.begin(),sortingArray.end() );
		double medianMag = sortingArray[ (N >> 1) ];
		// find the magnitude in the neighbourhoodMagnitudes
		unsigned int j = 0;
		while( neighbourhoodMagnitudes[j] != medianMag ) ++j;
		// set the data image pixel to the value of the temp image
		vtkIdType medianPoint = pointList->GetId( j );
		double *medianPixel = tempImage->GetPointData()->GetArray("Displacement")->GetTuple( medianPoint );
		this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i, medianPixel );
	}
}

/** A function to smooth the image using a weighted moving average using
 * a Gaussian kernel for weight calculation. The function is given R, 
 * the radius of points to consider; sigma, the standard deviation of the
 * of the Gaussian kernel; and mean, the mean of the Gaussian kernel in
 * terms of distance from the point being averaged (this will in 99.99%
 * of cases need to be set to 0).*/
void DisplacementWeightedMovingAverageFilter( double sigma, double mean )
{
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	
	// visit each point of the temp image, replace values in m_DataImage
	unsigned int nPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < nPoints; ++i ){
		
		double *newPixel = new double[3];
		newPixel = CalculateDisplacementWeightedMovingAverage( sigma, mean, i, tempImage );
		
		// set the data image pixel to the value of the new pixel
		this->m_DataImage->GetPointData()->GetArray("Displacement")->SetTuple( i, newPixel );
	}
}

/** A function to calculate the weighted average of the displacement 
 * stored in a pixel. This fucntion takes the standard deviation, sigma,
 * and the mean distance, mean, of the gaussian kernel that will be used
 * to caluculated the weights for the weighting function. It also take 
 * the mesh point ID, pointId and a pointer to the data image that 
 * should be used to caluclate the average. The point with ID pointID is
 * included in the calculation of the average. */ 
double* CalculateDisplacementWeightedMovingAverage( double sigma, double mean, unsigned int pointId, DataImagePointer image )
{
	if (sigma == 0) {
		return image->GetPointData()->GetArray("Displacement")->GetTuple( pointId );
	}
	
	vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New(); 
	this->GetNeighbouringPoints( pointId, pointList); // get connected points
	
	// Initialize variables for the loop
	double totalNumerator[3] = { 0 };
	double totalDenominator = 0;
	
	double *rPoint = new double[3];
	image->GetPoint( pointId, rPoint ); // reference point location	
	
	// visit each points in the list and use them to calculate the new point value
	for( unsigned int j = 0; j <= pointList->GetNumberOfIds(); ++j){
		vtkIdType cId = j == pointList->GetNumberOfIds() ? pointId : pointList->GetId( j ); // include the reference point at the very end
			
		double *cPoint = new double[3];
		image->GetPoint( cId, cPoint ); // current point location

		double distance = sqrt( pow(rPoint[0]-cPoint[0],2) + pow(rPoint[1]-cPoint[1],2) + pow(rPoint[2]-cPoint[2],2) );
		double weight = 1/(2.50662827*sigma)*pow(2.718281828,-pow((distance-mean),2)/(2*sigma*sigma)); // 1/(sqrt(2*pi)*sigma)*e^(-(X-mean)^2/(2*sigma^2)) (the Gaussian distribution)
		
		// add the weight to the total denominator
		totalDenominator = totalDenominator + weight;
		
		// add to the numerators
		double *tempPixel = new double[3];
		image->GetPointData()->GetArray("Displacement")->GetTuple( cId, tempPixel ); // the original data in the current point
		totalNumerator[0] = totalNumerator[0] + weight * tempPixel[0];
		totalNumerator[1] = totalNumerator[1] + weight * tempPixel[1];
		totalNumerator[2] = totalNumerator[2] + weight * tempPixel[2];
	}
	// use the total numerators and denominator to calculate the new pixel
	double *newPixel = new double[3];
	newPixel[0] = totalNumerator[0] / totalDenominator;
	newPixel[1] = totalNumerator[1] / totalDenominator;
	newPixel[2] = totalNumerator[2] / totalDenominator;
	
	return newPixel;
}

/** A function to smooth the image using a weighted moving average using
 * a Gaussian kernel for weight calculation. The function is given R, 
 * the radius of points to consider; sigma, the standard deviation of the
 * of the Gaussian kernel; and mean, the mean of the Gaussian kernel in
 * terms of distance from the point being averaged (this will in 99.99%
 * of cases need to be set to 0). Given a raduis of 0, the default value
 * of 3*sigma will be used*/
void StrainWeightedMovingAverageFilter( double sigma, double mean )
{
	// create a duplicate of the data image
	DataImagePointer tempImage = DataImagePointer::New();
	tempImage->DeepCopy(this->m_DataImage);
	
	// loop through the points of both the images
	unsigned int nPoints = tempImage->GetNumberOfPoints();
	for( unsigned int i = 0; i < nPoints; ++i ){
		
		double *newPixel = new double[9];
		newPixel = CalculateStrainWeightedMovingAverage( sigma, mean, i, tempImage );
		
		// set the data image pixel to the value of the new pixel
		this->m_DataImage->GetPointData()->GetArray("Strain")->SetTuple( i, newPixel );
	}
}

/** A function to calculate the weighted average of the strain 
 * stored in a pixel. This fucntion takes the standard deviation, sigma,
 * and the mean distance, mean, of the gaussian kernel that will be used
 * to caluculated the weights for the weighting function. It also take 
 * the mesh point ID, pointId and a pointer to the data image that 
 * should be used to caluclate the average. The point with ID pointID is
 * included in the calculation of the average. */
double* CalculateStrainWeightedMovingAverage( double sigma, double mean, unsigned int pointId, DataImagePointer image )
{
	if ( sigma == 0){
		return image->GetPointData()->GetArray("Stain")->GetTuple( pointId );
	}
	
	// get connected points
	vtkSmartPointer<vtkIdList> pointList = vtkSmartPointer<vtkIdList>::New(); // id list to hold ids of the closest points
	this->GetNeighbouringPoints( pointId, pointList );
	
	// Initialize vaiables for the loop
	double totalNumerator[9] = { 0 };
	double totalDenominator = 0;
	
	double *rPoint = new double[3]; // get reference point location
	image->GetPoint(pointId, rPoint );
	
	// visit each connected point in the list and use them to calculate the new point value
	for( unsigned int j = 0; j <= pointList->GetNumberOfIds(); ++j){
		vtkIdType cId = j == pointList->GetNumberOfIds() ? pointId : pointList->GetId( j ); // include the reference point at the very end

		double *cPoint = new double[3]; // get current point location
		image->GetPoint( cId, cPoint );
		
		double distance = sqrt( pow(*cPoint-*rPoint,2) + pow(*(cPoint+1)-*(rPoint+1),2) + pow(*(cPoint+2)-*(rPoint+2),2) );  //vector distance between the two
		double weight = 1/(2.50662827*sigma)*pow(2.718281828,-pow((distance-mean),2)/(2*sigma*sigma)); // 1/(sigma*sqrt(2*pi))*e^(-(X-mean)^2/(2*sigma^2)) (the Gaussian distribution)

		totalDenominator = totalDenominator + weight; // add the weight to the total denominator
		
		// add to the numerators
		double *tempPixel = image->GetPointData()->GetArray("Strain")->GetTuple( cId ); // the original data in the current point
		totalNumerator[0] = totalNumerator[0] + weight * tempPixel[0];
		totalNumerator[1] = totalNumerator[1] + weight * tempPixel[1];
		totalNumerator[2] = totalNumerator[2] + weight * tempPixel[2];
		totalNumerator[3] = totalNumerator[3] + weight * tempPixel[3];
		totalNumerator[4] = totalNumerator[4] + weight * tempPixel[4];
		totalNumerator[5] = totalNumerator[5] + weight * tempPixel[5];
		totalNumerator[6] = totalNumerator[6] + weight * tempPixel[6];
		totalNumerator[7] = totalNumerator[7] + weight * tempPixel[7];
		totalNumerator[8] = totalNumerator[8] + weight * tempPixel[8];
	}
	
	// use the total numerators and denominator to calculate the new pixel
	double *newPixel = new double[9];
	newPixel[0] = totalNumerator[0] / totalDenominator;
	newPixel[1] = totalNumerator[1] / totalDenominator;
	newPixel[2] = totalNumerator[2] / totalDenominator;
	newPixel[3] = totalNumerator[3] / totalDenominator;
	newPixel[4] = totalNumerator[4] / totalDenominator;
	newPixel[5] = totalNumerator[5] / totalDenominator;
	newPixel[6] = totalNumerator[6] / totalDenominator;
	newPixel[7] = totalNumerator[7] / totalDenominator;
	newPixel[8] = totalNumerator[8] / totalDenominator;
	
	return newPixel;
}

/** This function will return the region that the global registration 
 * registration should be conducted on. In the case of the DICMesh, this
 * region is the bounding box of the mesh. */
FixedImageRegionType GetGlobalRegistrationRegion()
{
	// restrict the Global registration to just the bounding box of the mesh
	double meshBBox[6];// = new double[6];
	this->m_DataImage->GetBounds( meshBBox );
	typename FixedImageType::PointType meshMinPt; // Get mesh minimum pt - this will be used for the registration region start index 
	meshMinPt[0] = meshBBox[0];
	meshMinPt[1] = meshBBox[2];
	meshMinPt[2] = meshBBox[4];
	double meshSize[3]; 				// Get mesh dimension - this will be used for the registration size
	meshSize[0] = meshBBox[1]-meshBBox[0];
	meshSize[1] = meshBBox[3]-meshBBox[2];
	meshSize[2] = meshBBox[5]-meshBBox[4];
		
	typename FixedImageType::IndexType fixedImageROIStart;
	fixedResampler->GetOutput()->TransformPhysicalPointToIndex(meshMinPt,fixedImageROIStart); // convert min point to start index
	
	typename FixedImageType::SpacingType fixedSpacing = fixedResampler->GetOutput()->GetSpacing(); // convert dimensinos to size in pixels
	typename FixedImageType::SizeType fixedImageROILengths;
	fixedImageROILengths[0] = (int)std::floor(meshSize[0]/fixedSpacing[0]);
	fixedImageROILengths[1] = (int)std::floor(meshSize[1]/fixedSpacing[1]);
	fixedImageROILengths[2] = (int)std::floor(meshSize[2]/fixedSpacing[2]);
	typename FixedImageType::RegionType fixedAnalysisRegion;
	fixedAnalysisRegion.SetIndex( fixedImageROIStart );
	fixedAnalysisRegion.SetSize( fixedImageROILengths );
	
	return fixedAnalysisRegion;
}

/** A function to set the downsample value used by the global
 * registration. The default value is 3. */
void SetGlobalRegistrationDownsampleValue( unsigned int value )
{
	if (m_GlobalRegDownsampleValue != value){
		this->m_GlobalRegDownsampleValue = value;
	}
}

/** A function to get the downsample value used by the global
 * registraion. */
unsigned int GetGlobalRegistrationDownsampleValue()
{
	return this->m_GlobalRegDownsampleValue;
}

/** A function to set the error tolerance in units of standard 
 * deviations from the average of the connected pixels. If a pixel is 
 * found to be farther from the mean than the value set, it will be
 * replaced by the ReplaceBadDisplacementPixels method. */
void SetDisplacementErrorTolerance(double tolerance)
{
	this->m_displacementErrorToleranceErrorTolerance = tolerance;
}

/** A function to get the error tolerance in units of standard 
 * deviations from the average of the connected pixels. If a pixel is 
 * found to be farther from the mean than the value set, it will be
 * replaced by the ReplaceBadDisplacementPixels method. */
double GetDisplacementErrorTollerance()
{
	return this->m_displacementErrorTolerance;
}

/** A function to set the error tolerance in units of standard 
 * deviations from the average of the connected pixels. If a pixel is 
 * found to be farther from the mean than the value set, it will be
 * replaced by the ReplaceBadStrainPixels method. */
void SetStrainErrorTolerance(double tolerance)
{
	this->m_strainErrorToleranceErrorTolerance = tolerance;
}

/** A function to get the error tolerance in units of standard 
 * deviations from the average of the connected pixels. If a pixel is 
 * found to be farther from the mean than the value set, it will be
 * replaced by the ReplaceBadStraintPixels method. */
double GetStrainErrorTolerance()
{
	return this->m_strainErrorTolerance;
}

/** A fucntion to set the error in image dimensions (eg. mm, inches). 
 * NOTE: Depriciated. This function may be removed at any time. Use at 
 * your OWN RISK.*/
void SetErrorRadius(double radius)
{
	this->m_errorRadius = radius;
}

/** A function to set the max metric value acceptable in registration.
 * NOTE: Depreciated. This function bay be removed at any time. Use at
 * your OWN RISK. */
void SetMaxMetricValue( double maxVal )
{
	if (m_maxMeticValue != maxVal){
		m_maxMeticValue = maxVal;
	}
}

/** A function to get the max metric value acceptable in registration.
 * NOTE: Depreciated. This function bay be removed at any time. Use at
 * your OWN RISK. */
double GetMaxMetricValue()
{
	return m_maxMeticValue;
}

/** A function that will return the points connected to the current
 * point, vtkIdType id, and returns them in vtkIdList idList. idList
 * is reset to empty by this method before being filled.*/
void GetNeighbouringPoints( vtkIdType id, vtkSmartPointer<vtkIdList> idList )
{
	idList->Reset(); // reset the ID list to avoid mistakes
	
	// Get the cell IDs that 'id' is part of
	vtkSmartPointer<vtkIdList> cCellList = vtkSmartPointer<vtkIdList>::New();
	this->m_DataImage->GetPointCells( id, cCellList );
	
	// visit each cell
	for( unsigned int i = 0; i < cCellList->GetNumberOfIds(); ++i ){
		
		// visit each line in each cell
		for ( int j = 0; j < this->m_DataImage->GetCell( cCellList->GetId( i ) )->GetNumberOfEdges(); ++j){
			
			// if 'id' is the first in the line, save the second id as a unique id
			if( this->m_DataImage->GetCell(cCellList->GetId(i))->GetEdge(j)->GetPointId(0) == id ){
				idList->InsertUniqueId( this->m_DataImage->GetCell(cCellList->GetId(i))->GetEdge(j)->GetPointId(1) );
				continue;
			}
			// if 'id' is the second in the line, save the first id as a unique id
			if( this->m_DataImage->GetCell(cCellList->GetId(i))->GetEdge(j)->GetPointId(1) == id ){
				idList->InsertUniqueId( this->m_DataImage->GetCell(cCellList->GetId(i))->GetEdge(j)->GetPointId(0) );
				continue;
			}
		}
	}
}

private:

DataImagePointer			m_DataImage;
vtkSmartPointer<vtkPKdTree>	m_KDTree;
double						m_errorRadius;
double						m_displacementErrorTolerance;
double						m_strainErrorTolerance;
double						m_maxMeticValue; // this because I'm using the normalized x-correlation coefficient metric
vtkSmartPointer<vtkIdList>	m_pointsList;
RegistrationParametersType	m_GlobalRegistrationParameters;
unsigned int				m_GlobalRegDownsampleValue;
	
}; // end class DICMesh


#endif // DICMESH_H

