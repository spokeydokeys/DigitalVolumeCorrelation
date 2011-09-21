//      DIC.cxx
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

#ifndef DIC_H
#define DIC_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "itkImage.h"
#include "itkVector.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// includes for the registration
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkNormalizeImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkTranslationTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCenteredEuler3DTransform.h"
#include <itkCenteredTransformInitializer.h>


template <typename TFixedImage, typename TMovingImage>
class DIC
{
public:
/** Types **/

/** Type of the Fixed Image */
typedef				TFixedImage						FixedImageType;
typedef	typename	FixedImageType::ConstPointer	FixedImageConstPointer;
typedef	typename	FixedImageType::RegionType		FixedImageRegionType;
typedef	typename	FixedImageType::Pointer			FixedImagePointer;

/** Type of the Moving Image */
typedef				TMovingImage					MovingImageType;
typedef	typename	MovingImageType::ConstPointer	MovingImageConstPointer;
typedef	typename	MovingImageType::RegionType		MovingImageRegionType;
typedef	typename	MovingImageType::Pointer		MovingImagePointer;

/** Type to hold the number of regions. */
typedef	unsigned int				NumberOfInterrogationRegionsType;

/** Type of the scaling factor for the Moving Image IRRadius. For
 * efficiency only a section of the moving image is considered for
 * search.  If "UseWholeMovingImage" is set to true, the whole moving
 * image will be held in memory. */
typedef double						FixedImageIRMultiplierType;

bool								UseWholeMovingImage;  // not yet implemented

/** Type of the Fixed Region List. */
typedef	std::vector< FixedImageRegionType* >									FixedImageRegionListType;

/** Type of the Moving Region List. */
typedef	std::vector< MovingImageRegionType* >									MovingImageRegionListType;

/** Type of the Fixed Image ROI filter. */
typedef	itk::RegionOfInterestImageFilter< FixedImageType, FixedImageType >		FixedROIFilterType;
typedef	typename	FixedROIFilterType::Pointer									FixedROIFilterPointer;
typedef itk::ImageFileReader< FixedImageType >									FixedImageReaderType;
typedef typename	FixedImageReaderType::Pointer								FixedImageReaderPointer;

/** Type of the Moving Image ROI filter. */
typedef	itk::RegionOfInterestImageFilter< MovingImageType, MovingImageType >	MovingROIFilterType;
typedef	typename	MovingROIFilterType::Pointer								MovingROIFilterPointer;
typedef itk::ImageFileReader< MovingImageType >									MovingImageReaderType;
typedef typename	MovingImageReaderType::Pointer								MovingImageReaderPointer;

/** Type of the Image Registration Method */
typedef	itk::ImageRegistrationMethod< FixedImageType, MovingImageType>			ImageRegistrationMethodType;
typedef	typename	ImageRegistrationMethodType::Pointer						ImageRegistrationMethodPointer;

typedef itk::NormalizedCorrelationImageToImageMetric< FixedImageType, MovingImageType >	MetricType;
typedef typename	MetricType::Pointer											MetricTypePointer;

typedef itk::RegularStepGradientDescentOptimizer								OptimizerType;
typedef typename	OptimizerType::Pointer										OptimizerTypePointer;

typedef itk::CenteredEuler3DTransform< double >									TransformType;
typedef	typename	TransformType::Pointer										TransformTypePointer;

typedef itk::LinearInterpolateImageFunction< MovingImageType, double >			InterpolatorType;
typedef typename	InterpolatorType::Pointer									InterpolatorTypePointer;

typedef itk::CenteredTransformInitializer<TransformType,FixedImageType,MovingImageType>	TransformInitializerType;
typedef typename	TransformInitializerType::Pointer							TransformInitializerTypePointer;

/** Methods **/
/** Constructor **/
DIC()
{
	m_FixedImage			= 0; // must be provided by user
	m_MovingImage			= 0; // must be provided by user
	m_IRRadius				= 0; // must be provided by user
	//~ m_Registration			= 0; // mut be provided by user
	
	
	m_CurrentFixedImage		= 0;
	m_CurrentMovingImage	= 0;

	UseWholeMovingImage		= false;
	m_FixedIRMult			= 2.5; // default size of the moving IR is 2.5 times the size of the fixed IR
	
	m_FixedROIFilter		= FixedROIFilterType::New();
	m_MovingROIFilter		= MovingROIFilterType::New();
	
	// setup the registration
	m_Metric				= MetricType::New();	
	m_Optimizer				= OptimizerType::New();
	m_Transform				= TransformType::New();
	m_Interpolator			= InterpolatorType::New();
	m_Registration			= ImageRegistrationMethodType::New();
	m_Registration->SetTransform( m_Transform );
	m_Registration->SetInterpolator( m_Interpolator );
	m_Registration->SetOptimizer( m_Optimizer );
	m_Registration->SetMetric( m_Metric );
	m_TransformInitializer	= TransformInitializerType::New();

}

/** Destructor **/
~DIC() {}

/** Set the Fixed image. */
void SetFixedImage( FixedImageConstPointer fixedImage )
{
	if (this->m_FixedImage.GetPointer() != fixedImage )
	{
		this->m_FixedImage = fixedImage;
	}
}

/** Get the Fixed Image. */
FixedImageConstPointer GetFixedImage()
{
	return this->m_FixedImage;
}

/** Set the Moving image. */
void SetMovingImage( MovingImageConstPointer movingImage)
{
	if (this->m_MovingImage.GetPointer() != movingImage )
	{
		this->m_MovingImage = movingImage;
	}	
}

/** Get the Moving Image. **/
MovingImageConstPointer GetMovingImage()
{
	return this->m_MovingImage;
}

/** Set the Interrogation Region Radius. */
void SetInterrogationRegionRadius( unsigned int radius )
{
	if (this->m_IRRadius != radius )
	{
		this->m_IRRadius = radius;
	}
}

/** Get the current Interrogation Region Radius. */
unsigned int GetInterrogationRegionRadius()
{
	return this->m_IRRadius;
}

/** Set the Registration Method. */
void SetRegistrationMethod( ImageRegistrationMethodPointer registrationMethod )
{
	if (this->m_Registration.GetPointer() != registrationMethod)
	{
		this->m_Registration = registrationMethod;
	}
}

/** Get the Registration method. */
ImageRegistrationMethodPointer GetRegistrationMethod()
{
	return this->m_Registration;
}

/** A function to return the i'th fixed image region. */
FixedImageRegionType* GetFixedImageRegionFromIndex( unsigned int regionNumber)
{

	if (this->m_FixedImageRegionList.empty())
	{
		std::stringstream msg("");
		msg << "The fixed region list is currently empty."<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	
	return this->m_FixedImageRegionList.at( regionNumber );
}

/** a function to get the fixed image region list. Returns the list pointer. **/
FixedImageRegionListType* GetFixedImageRegionList()
{
	return &this->m_FixedImageRegionList;
}

/** A function to return the i'th moving image region. */
MovingImageRegionType* GetMovingImageRegionFromIndex( unsigned int regionNumber )
{
	if (this->m_MovingImageRegionList.empty())
	{
		std::stringstream msg("");
		msg << "The moving region list is currently empty."<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
		
	return this->m_MovingImageRegionList.at( regionNumber );
}

/** a function to add a region to the end of the fixed image region list. **/
void PushRegionOntoFixedImageRegionList( FixedImageRegionType *region )
{
	this->m_FixedImageRegionList.push_back( region );

}

/** This function takes the real-world location and returns the region
 * of the moving image based on the IRRadius. */
void GetMovingImageRegionFromLocation( MovingImageRegionType *region, double *centerLocation )
{
	typename MovingImageType::IndexType						centerIndex;
	typename MovingImageType::IndexType						startIndex;
	typename MovingImageRegionType::SizeType				regionSize;
	typename MovingImageRegionType::SizeType::SizeValueType	regionDiameter;	

	this->m_MovingImage->TransformPhysicalPointToIndex( (typename MovingImageType::PointType) centerLocation, centerIndex );
	regionDiameter = this->m_IRRadius*2+1;
	regionSize.Fill(regionDiameter);
	startIndex[0] = centerIndex[0] - this->m_IRRadius;
	startIndex[1] = centerIndex[1] - this->m_IRRadius;
	startIndex[2] = centerIndex[2] - this->m_IRRadius;
	
	region->SetSize(regionSize); // set the parameters of the region
	region->SetIndex(startIndex);
	//std::cout<<"Moving region: "<<*region<<std::endl;
	if ( !this->IsMovingRegionValid( region) ) // if the region is invalid, change the parameters
	{
		
		this->FixImageRegion( region, m_MovingImage );
		//std::cout<<"Fixed moving region: "<<*region<<std::endl;
	}
}
 
/** This function takes the real-world location and returns the region
 * of the fixed image based on the IRRadus. */
void GetFixedImageRegionFromLocation( FixedImageRegionType *region, double *centerLocation )
{
	typename FixedImageType::IndexType						centerIndex;
	typename FixedImageType::IndexType						startIndex;
	typename FixedImageRegionType::SizeType					regionSize;
	typename FixedImageRegionType::SizeType::SizeValueType	regionDiameter;	

	this->m_FixedImage->TransformPhysicalPointToIndex( (typename FixedImageType::PointType) centerLocation, centerIndex );
	regionDiameter = this->m_IRRadius*this->m_FixedIRMult*2+1;
	regionSize.Fill(regionDiameter);
	startIndex[0] = centerIndex[0] - this->m_IRRadius*this->m_FixedIRMult;
	startIndex[1] = centerIndex[1] - this->m_IRRadius*this->m_FixedIRMult;
	startIndex[2] = centerIndex[2] - this->m_IRRadius*this->m_FixedIRMult;
	
	region->SetSize(regionSize); // set the parameters of the region
	region->SetIndex(startIndex);
	//std::cout<<"Fixed region: "<<*region<<std::endl;
	if ( !this->IsFixedRegionValid( region ) ) // if the region is invalid, change the parameters
	{
		
		this->FixImageRegion( region, m_FixedImage );
		//std::cout<<"Fixed fixed region: "<<*region<<std::endl;
	}
}

/** This funciton will determine if a given fixed image region is valid.*/
bool IsFixedRegionValid( FixedImageRegionType *region )
{
	//~ return this->m_FixedImage->GetLargestPossibleRegion().IsInside( *region );
	typename FixedImageRegionType::SizeType	size = region->GetSize();
	if( !this->m_FixedImage->GetLargestPossibleRegion().IsInside( *region ) || size[0] < 4 || size[1] < 4 || size[2] < 4){
		return false;
	}
	
	return true;
}

/** This function will determine if a given moving image region is valid. */
bool IsMovingRegionValid( MovingImageRegionType *region )
{
	//~ return this->m_MovingImage->GetLargestPossibleRegion().IsInside( *region );
	typename FixedImageRegionType::SizeType	size = region->GetSize();
	if( !this->m_MovingImage->GetLargestPossibleRegion().IsInside( *region ) || size[0] < 4 || size[1] < 4 || size[2] < 4){
		return false;
	}
	
	return true;	
}

//~ /** This function will modify a fixed image region to make it valid.*/
//~ void FixFixedImageRegion( FixedImageRegionType *region )
//~ {
	//~ FixedImageRegionType						fixedImageRegion		= this->m_FixedImage->GetLargestPossibleRegion();	
	//~ typename FixedImageRegionType::IndexType	fixedImageRegionIndex 	= fixedImageRegion.GetIndex();
	//~ typename FixedImageRegionType::SizeType		fixedImageRegionSize 	= fixedImageRegion.GetSize();
	//~ unsigned int								fixedImageDimension 	= this->m_FixedImage->GetImageDimension();
	//~ 
	//~ typename FixedImageRegionType::SizeType		regionSize				= region->GetSize();
	//~ typename FixedImageRegionType::IndexType	regionIndex				= region->GetIndex();
	//~ 
	//~ for (unsigned int i = 0; i<fixedImageDimension; ++i)
	//~ {
		//~ if ( regionIndex[i] < fixedImageRegionIndex[i] ) // if region starts before the image
		//~ {
			//~ regionSize[i] = regionSize[i] - (fixedImageRegionIndex[i] - regionIndex[i]); // change size to that inside the image
			//~ regionIndex[i] = fixedImageRegionIndex[i]; // change the index to the start of the image
			//~ regionSize[i] > m_FixedIRMult*m_IRRadius ?  : regionSize[i] = m_FixedIRMult*m_IRRadius; // the gaussian filter used in interpolation needs >= 4 in each direction.
		//~ }
		//~ if ( regionIndex[i]  > (int)(fixedImageRegionIndex[i] + fixedImageRegionSize[i]) ) // if the region starts after the image finishes, throw a warning and use the closest part of the image.
		//~ {
			//~ std::stringstream msg("");
			//~ msg<<"Warning: The fixed image region starts after the end of the moving image.  Using closest valid region."<<std::endl<<"Index: ["<<regionIndex[0]<<", "<<regionIndex[1]<<", "<<regionIndex[2]<<"]"<<std::endl<<"Size: ["<<regionSize[0]<<", "<<regionSize[1]<<", "<<regionSize[2]<<"]"<<std::endl;
			//~ regionIndex[i] = fixedImageRegionSize[i] - (m_FixedIRMult*m_IRRadius+1);
			//~ regionSize[i] = m_FixedIRMult*m_IRRadius;
		//~ }
		//~ if ( (regionIndex[i] + regionSize[i]) > (fixedImageRegionIndex[i] + fixedImageRegionSize[i]) ) // if the region finishes after the image
		//~ {
			//~ regionSize[i] = (fixedImageRegionIndex[i] + fixedImageRegionSize[i]) - (regionIndex[i]+1); // change size to that inside the image
			//~ if (regionSize[i] < m_FixedIRMult*m_IRRadius){
				//~ regionIndex[i] = fixedImageRegionSize[i] - (m_FixedIRMult*m_IRRadius+1);
				//~ regionSize[i] = m_FixedIRMult*m_IRRadius;
			//~ } // the gaussian filter used in interpolation needs >= 4 in each direction, index is 0 referenced, so needs -5
		//~ }
		//~ 
	//~ }
	//~ 
	//~ region->SetSize( regionSize );
	//~ region->SetIndex( regionIndex );
//~ }
//~ 
//~ /** This function will modify a moving image region to make it valid.*/
//~ void FixMovingImageRegion( MovingImageRegionType *region )
//~ {
	//~ MovingImageRegionType						movingImageRegion		= this->m_MovingImage->GetLargestPossibleRegion();
	//~ typename MovingImageRegionType::IndexType	movingImageRegionIndex	= movingImageRegion.GetIndex();
	//~ typename MovingImageRegionType::SizeType	movingImageRegionSize	= movingImageRegion.GetSize();
	//~ unsigned int								movingImageDimension	= this->m_MovingImage->GetImageDimension();
	//~ 
	//~ typename MovingImageRegionType::SizeType	regionSize 				= region->GetSize();
	//~ typename MovingImageRegionType::IndexType	regionIndex				= region->GetIndex();
	//~ 
	//~ for (unsigned int i = 0; i<movingImageDimension; ++i)
	//~ {
		//~ if ( regionIndex[i] < movingImageRegionIndex[i] ) // if region starts before the image
		//~ {
			//~ regionSize[i] = regionSize[i]-(movingImageRegionIndex[i]-regionIndex[i])>m_IRRadius ? regionSize[i]-(movingImageRegionIndex[i]-regionIndex[i]) : m_IRRadius; // change size to that inside the image, but use at least the perscribed IRRadius
			//~ regionIndex[i] = movingImageRegionIndex[i]; // change the index to the start of the image
			//~ regionSize[i] = regionSize[i] > m_IRRadius ? : regionSize[i] = m_IRRadius; // the gaussian filter used in interpolation needs >= 4 in each direction.
		//~ }
		//~ if ( regionIndex[i] > (int)(movingImageRegionIndex[i] + movingImageRegionSize[i]) ) // if the region starts after the image finishes, throw a warning and use the closest part of the image.
		//~ {
			//~ std::stringstream msg("");
			//~ msg<<"Warning: The moving image region starts after the end of the moving image.  Using closest valid region."<<std::endl<<"Index: ["<<regionIndex[0]<<", "<<regionIndex[1]<<", "<<regionIndex[2]<<"]"<<std::endl<<"Size: ["<<regionSize[0]<<", "<<regionSize[1]<<", "<<regionSize[2]<<"]"<<std::endl;
			//~ regionIndex[i] = movingImageRegionSize[i] - (m_IRRadius+1);
			//~ regionSize[i] = m_IRRadius;
		//~ }		
		//~ if ( (regionIndex[i] + regionSize[i]) > (movingImageRegionIndex[i] + movingImageRegionSize[i]) ) // if the region finishes after the image
		//~ {
			//~ regionSize[i] = (movingImageRegionIndex[i] + movingImageRegionSize[i]) - (regionIndex[i]+1); // change size to that inside the image
			//~ if (regionSize[i] < m_IRRadius){
				//~ regionIndex[i] = movingImageRegionSize[i] - (m_IRRadius+1);
				//~ regionSize[i] = m_IRRadius;
			//~ } // the gaussian filter used in interpolation needs >= 4 in each direction, index is 0 referenced, so needs -5
		//~ }
//~ 
	//~ }
	//~ region->SetSize( regionSize );
	//~ region->SetIndex( regionIndex );
//~ }

/** This function will modify a moving image region to make it valid.*/
void FixImageRegion( MovingImageRegionType *region, MovingImageConstPointer image )
{
	MovingImageRegionType						movingImageRegion		= image->GetLargestPossibleRegion();
	typename MovingImageRegionType::IndexType	movingImageRegionIndex	= movingImageRegion.GetIndex();
	typename MovingImageRegionType::SizeType	movingImageRegionSize	= movingImageRegion.GetSize();
	unsigned int								movingImageDimension	= this->m_MovingImage->GetImageDimension();
	
	typename MovingImageRegionType::SizeType	regionSize 				= region->GetSize();
	typename MovingImageRegionType::IndexType	regionIndex				= region->GetIndex();
	
	for (unsigned int i = 0; i<movingImageDimension; ++i)
	{
		if ( regionIndex[i] < movingImageRegionIndex[i] ) // if region starts before the image
		{
			regionSize[i] = regionSize[i]-(movingImageRegionIndex[i]-regionIndex[i])>m_IRRadius ? regionSize[i]-(movingImageRegionIndex[i]-regionIndex[i]) : m_IRRadius; // change size to that inside the image, but use at least the perscribed IRRadius
			regionIndex[i] = movingImageRegionIndex[i]; // change the index to the start of the image
		}
		if ( regionIndex[i] > (int)(movingImageRegionIndex[i] + movingImageRegionSize[i]) ) // if the region starts after the image finishes, throw a warning and use the closest part of the image.
		{
			std::stringstream msg("");
			msg<<"Warning: The moving image region starts after the end of the moving image.  Using closest valid region."<<std::endl<<"Index: ["<<regionIndex[0]<<", "<<regionIndex[1]<<", "<<regionIndex[2]<<"]"<<std::endl<<"Size: ["<<regionSize[0]<<", "<<regionSize[1]<<", "<<regionSize[2]<<"]"<<std::endl;
			regionIndex[i] = movingImageRegionSize[i] - (m_IRRadius+1);
			regionSize[i] = m_IRRadius;
		}		
		if ( (regionIndex[i] + regionSize[i]) > (movingImageRegionIndex[i] + movingImageRegionSize[i]) ) // if the region finishes after the image
		{
			regionSize[i] = (movingImageRegionIndex[i] + movingImageRegionSize[i]) - (regionIndex[i]+1); // change size to that inside the image
			if (regionSize[i] < m_IRRadius){
				regionIndex[i] = movingImageRegionSize[i] - (m_IRRadius+1);
				regionSize[i] = m_IRRadius;
			} // use a minimum size of the IRRadius
		}

	}
	region->SetSize( regionSize );
	region->SetIndex( regionIndex );
}

/** A function to get the moving image region list. Returns the list pointer **/
MovingImageRegionListType* GetMovingImageRegionList()
{
	return &this->m_MovingImageRegionList;
}

/** A function to add a region to the end of the moving image region list. **/
void PushRegionOntoMovingImageRegionList( MovingImageRegionType *region )
{
	this->m_MovingImageRegionList.push_back( region );
}

/** A function to get the Fixed ROI as an image from the ROI filter.*/
FixedImagePointer GetFixedROIAsImage( FixedImageRegionType *desiredRegion )
{
	this->m_FixedROIFilter->SetInput( this->m_FixedImage );
	this->m_FixedROIFilter->SetRegionOfInterest( *desiredRegion );
	
	try
	{
		this->m_FixedROIFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::stringstream msg("");
		msg <<"Exception caught updating GetFixedROIAsImage method" << std::endl << "Exception message" << 
			std::endl << err << std::endl << std::endl;
		this->WriteToLogfile( msg.str() );		
		std::abort();
	}
	return this->m_FixedROIFilter->GetOutput();
}

/** A function to get the Moving ROI as an image from the ROI filter.*/
MovingImagePointer GetMovingROIAsImage( MovingImageRegionType *desiredRegion )
{
	this->m_MovingROIFilter->SetInput( this->m_MovingImage );
	this->m_MovingROIFilter->SetRegionOfInterest( *desiredRegion );
	
	try
	{
		this->m_MovingROIFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::stringstream msg("");
		msg << "Exception caught updating GetMovingROIAsImage method" << std::endl << "Exception message" <<
			std::endl<<err<<std::endl<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	
	return this->m_MovingROIFilter->GetOutput();
}

/** A function to set the fixed roi image as the registration image.*/
void SetFixedROIImage( FixedImagePointer roiImage )
{
	if (this->m_CurrentFixedImage.GetPointer() != roiImage.GetPointer() )
	{
		this->m_CurrentFixedImage = roiImage;
	}
}

/** A function to set the moving roi image as the registration image.*/
void SetMovingROIImage( MovingImagePointer roiImage )
{
	if (this->m_CurrentMovingImage.GetPointer() != roiImage.GetPointer() )
	{
		this->m_CurrentMovingImage = roiImage;
	}
}

/** A function to get the registration of the two regions. */
void UpdateRegionRegistration()
{
	this->m_Registration->SetFixedImage( this->m_CurrentFixedImage );
	this->m_Registration->SetMovingImage( this->m_CurrentMovingImage );

	try
	{
		std::stringstream msg("");
		msg << "Updating Registration!"<<std::endl;
		this->WriteToLogfile( msg.str() );
		this->m_Registration->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::stringstream msg("");
		msg <<"Exception caught updating GetRegionRegistration method"<<std::endl<<"Exception message"<<
			std::endl<<err<<std::endl<<std::endl;
		this->WriteToLogfile( msg.str() );
				
		typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::OptimizerType::Pointer optimizer;
		optimizer = this->m_Registration->GetOptimizer();
		msg.str(" ");
		msg << optimizer<<std::endl;
		this->WriteToLogfile( msg.str() );
		
		typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::FixedImageType::ConstPointer	fixedImage;
		fixedImage = this->m_Registration->GetFixedImage();
		msg.str(" ");
		msg << fixedImage << std::endl;
		this->WriteToLogfile( msg.str() );	
				
		typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::MovingImageType::ConstPointer	movingImage;
		movingImage = this->m_Registration->GetMovingImage();
		msg.str(" ");
		msg << movingImage << std::endl;
		this->WriteToLogfile( msg.str() );
	}
}

/** A fucntion that modifies an array of three doubles to containe
 * the final displacement from the registration.
 * This function assumes that the last three parameters of the transform
 * parameters are the displacement. Depending on your chosen transfrom,
 * this may not be the case and you should check the ITK documentation
 * to see.
 * TODO: make a vector of equal length to the parameters vector called
 * transformDisplacements, with type bool.  1's define locations of
 * displacement data in the parameters.**/
void GetLastDisplacement( double *pixelData )
{
	typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType finalParameters = this->m_Registration->GetLastTransformParameters();
	unsigned int nParameters = this->m_Registration->GetTransform()->GetNumberOfParameters();
	*pixelData 		= finalParameters[nParameters-3];
	*(pixelData +1)	= finalParameters[nParameters-2];
	*(pixelData +2)	= finalParameters[nParameters-1];
}

/** A function that modifies a single double to contain the optmizer
 * value from the last optimization.**/
void GetLastOptimizer( double *optData )
{
	typename DIC<FixedImageType, MovingImageType>::ImageRegistrationMethodType::ParametersType finalParameters = this->m_Registration->GetLastTransformParameters();
	typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::OptimizerType::Pointer optimizer = this->m_Registration->GetOptimizer();
	*optData = optimizer->GetValue( finalParameters );
}

/** A function that sets the initial displacement in the registration. 
 * This function assumes that the last three parameters of the
 * transform parameters are the displacement. Depending on your chosen
 * transform, this may not be the case and you should check the itk 
 * documentation to see. 
 * TODO: make a vector of equal length to the parameters vector called
 * transformDisplacements, with type bool.  1's define locations of
 * displacement data in the parameters.**/
void SetInitialDisplacement(double *pixelData)
{
	typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType initialParameters;
	initialParameters = this->m_Registration->GetInitialTransformParameters();
	initialParameters[6] = *pixelData;
	initialParameters[7] = *(pixelData +1);
	initialParameters[8] = *(pixelData +2);
	
	this->m_Registration->SetInitialTransformParameters( initialParameters );
}

/** A function to set the rotation centre of the registration. 
 * This function assumes that the 4th 5th and 6th parameters of the
 * transform parameters are the rotation centre. Depending on your chosen
 * transform, this may not be the case and you should check the itk 
 * documentation to see.  
 * TODO: make a vector of equal length to the parameters vector called
 * transformRotationCenters, with type bool.  1's define locations of
 * rotation centers in the parameters. **/
void SetRotationCenter( double *regionCenter )
{
	typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType initialParameters;
	initialParameters = this->m_Registration->GetInitialTransformParameters();
	initialParameters[3] = *regionCenter;
	initialParameters[4] = *(regionCenter+1);
	initialParameters[5] = *(regionCenter+2);
	
	this->m_Registration->SetInitialTransformParameters( initialParameters );
}

/** A function to set the transform to the identity transform. 
 * This function currently is only for the rotational transform.
 * TODO: make local a variable called identityTransform that must be
 * set by the user and will be generic.*/
void SetTransformToIdentity()
{
	typename DIC<FixedImageType,MovingImageType>::ImageRegistrationMethodType::ParametersType initialParameters;
	initialParameters = this->m_Registration->GetInitialTransformParameters();
	initialParameters[0] = 0;
	initialParameters[1] = 0;
	initialParameters[2] = 0;
	initialParameters[3] = 0;
	initialParameters[4] = 0;
	initialParameters[5] = 0;
	initialParameters[6] = 0;
	initialParameters[7] = 0;
	initialParameters[8] = 0;

	this->m_Registration->SetInitialTransformParameters( initialParameters );
}

/** A function to set the logfile. */
void SetLogfileName( std::string logfile )
{
	this->m_LogfileName = logfile;	
}

/** A function that returns the logfile name as a stringn */
std::string GetLogfileName()
{
	return this->m_LogfileName;
}

/** A function to write string data to a log file. */
void WriteToLogfile( std::string characters)
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

/** A function to set the output directory .*/
void SetOuputDirectory( std::string directory )
{
	this->m_OutputDirectory = directory;
}

/** A function to get the output directory.*/
std::string GetOutputDirectory()
{
	return this->m_OutputDirectory;
}

/** A function to set the transform. */
void SetTranform( TransformTypePointer transform )
{
	if (m_Transform.GetPointer() != transform){
		this->m_Transform = transform;
	}
}

/** A function to get the transform. */
TransformTypePointer GetTransform()
{
	return this->m_Transform;
}

/** A method to set the optimizer. */
void SetOptimizer( OptimizerTypePointer optimizer )
{
	if (m_Optimizer.GetPointer() != optimizer ){
		this->m_Optimizer = optimizer;
	}
}

/** A method to get the optimizer. */
OptimizerTypePointer GetOptimizer()
{
	return this->m_Optimizer;
}
	
protected:

FixedImageConstPointer				m_FixedImage;
MovingImageConstPointer				m_MovingImage;

FixedImagePointer					m_CurrentFixedImage;
MovingImagePointer					m_CurrentMovingImage;

FixedImageRegionListType			m_FixedImageRegionList;
MovingImageRegionListType			m_MovingImageRegionList;

FixedROIFilterPointer				m_FixedROIFilter;
MovingROIFilterPointer				m_MovingROIFilter;

unsigned int						m_IRRadius;
NumberOfInterrogationRegionsType	m_NoIRs;
FixedImageIRMultiplierType			m_FixedIRMult;

ImageRegistrationMethodPointer		m_Registration;
MetricTypePointer					m_Metric;
OptimizerTypePointer				m_Optimizer;
TransformTypePointer				m_Transform;
InterpolatorTypePointer				m_Interpolator;
TransformInitializerTypePointer		m_TransformInitializer;

std::string							m_LogfileName;
std::string							m_OutputDirectory;
}; // end class DIC

#endif // DIC_H
