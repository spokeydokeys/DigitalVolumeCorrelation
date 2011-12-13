//      DICGrid.cxx
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

#ifndef DICGRID_H
#define DICGRID_H

#include <cstring>
#include <ctime>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkStructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellDerivatives.h>
#include <vtkCellDataToPointData.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageData.h>

#include "itkImageFileReader.h"
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMap.h>
#include "DIC.cxx"



template<typename TFixedImage, typename TMovingImage, typename TMaskImage>
class DICGrid : public DIC<TFixedImage, TMovingImage>
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
typedef typename DIC<TFixedImage,TMovingImage>::TransformInitializerType	TransformInitializerType;
typedef typename DIC<TFixedImage,TMovingImage>::TransformType				TransformType;

typedef			 TMaskImage							MaskImageType;
typedef	typename MaskImageType::Pointer				MaskImagePointer;
typedef	typename MaskImageType::ConstPointer		MaskImageConstPointer;

typedef		vtkSmartPointer<vtkImageData>			DataImagePointer;
typedef		vtkSmartPointer<vtkPoints>				DataImagePointsPointer;
typedef		vtkSmartPointer<vtkDoubleArray>			DataImagePixelPointer;

typedef		double									RegionOverlapType;
typedef		double									MaskRatioType;

typedef	itk::ShapeLabelObject< unsigned long, FixedImageType::GetImageDimension >	LabelObjectType;
typedef	typename LabelObjectType::Pointer											LabelObjectPointer;
typedef typename itk::LabelMap< LabelObjectType >									LabelMapType;
typedef typename LabelMapType::Pointer												LabelMapPointer;


/* Methods. **/

/** Constructor **/
DICGrid()
{
	m_DataImage = 0; // must be provided by user or read in using the ReadMeshFromGmshFile method
	m_errorRadius = 4; // initialize the error search radius to 3 image units
	m_displacementErrorTolerance = 2; // difference of a pixel from its neighbours to be considered erronious, in standard deviations from the mean
	m_strainErrorTolerance = 1;	
	m_pointsList = vtkSmartPointer<vtkIdList>::New(); // the points list for analysis
	m_maxMeticValue = -0.00; // TODO: make this setable using a method
	m_maskImage = 0; // Must be provided by the user
	m_regionOverlap = 0.5; // default region overlap
	m_maskRatio = 0.3; // the minimum value of region/mask overlap to evaluate the region
	m_labelMap = 0; // to be calculated from the mask image
}

void CreateEmptyDataImage()
{
	// get the centre of the data image
	typename LabelObjectType::CentroidType	center = this->GetMaskImageLabel()->GetCentroid();
	typename LabelMapType::IndexType centerIndex;
	this->GetLabelMap()->TransformPhysicalPointToIndex( center, centerIndex );
	
	// get the point spacing in the image
	double spacing = (this->GetInterrogationRegionRadius()*2+1)*this->GetRegionOverlap();
	
	// get the number of points to the edge of the label mask
	typename MaskImageType::RegionType	labelMapRegion = this->GetLabelMap()->GetRegion();
	 
	DataImagePointer	dataImage = DataImagePointer::New();
	
	
	
	
	
	
	// Find the middle of the image.
	// use itkImageMomentsCalculator to calculate the mask center
	// example in: /media/data/ITK-Source/InsightToolkit-3.20.0/Code/Common/itkCenteredTransformInitializer.txx
	
	// Find the number of Regions in each direction in the image
	//		use m_regionOverlap and m_IRRadius to define an image spacing.
	// create vtkStructuredGrid to with the correct size, spacing and origin
}

/** Set the region overlap. For practicality, the region overlap is limited
 * to values < 0.95. If the region overlap is >= 0.95, or if the value
 * is not set, this method returns false. If the method is set or if the
 * value was not changed due to similarity this method returns true. */
bool SetRegionOverlap( RegionOverlapType overlap )
{
	if (overlap >= 0.95 ){
		std::stringstream msg("");
		msg<<"Region overlap must be < 0.95. Please check overlap and try again."<<std::endl;
		this->WriteToLogfile( msg.str() );
		return false;
	}
	if (this->m_regionOverlap != overlap){
		this->m_regionOverlap = overlap;
		return true;
	}
	else if( this->m_regionOverlap == overlap){
		return true;
	}
	else{
		return false;
	}
}

RegionOverlapType GetRegionOverlap()
{
	return this->m_regionOverlap;
}

void CalculateInitialMovingImageRgionList()
{
	// loop through the points in the data image
	// check if the region overlaps the mask image
	//		use m_maskRatio
	// if it does not overlap the mask use vtkStructuredGrid::BlankPoint to turn off the point
	// it is does overlap, add the region, defined by m_IRRadius to the region list
	// to speed things up, perform the following psudo code:
	//  if(center of region is contained in object => store region as valid
	//  if(!one of the corners is in object) => skip as invalid region
	//  if(!ratio of region in image > m_maskRatio) => skip as invalid region
}

void CalculateInitialFixedImageRegionList()
{
	// loop through the points in the data image
	// check if the region overlaps the mask image
	//		use m_maskRatio
	// if it does not overlap the mask use vtkStructuredGrid::BlankPoint to turn off the point
	// it is does overlap, add the region, defined by m_IRRadius to the region list
	// to speed things up, perform the following psudo code:
	//  if(center of region is contained in object => store region as valid
	//  if(!one of the corners is in object) => skip as invalid region
	//  if(!ratio of region in image > m_maskRatio) => skip as invalid region
}

FixedImageRegionType GetGlobalRegistrationRegion()
{
	// return the bounding region of the largest labeled region in the mask image
}

void SetMaskRatio(double ratio)
{
	if( ratio != m_maskRatio ){
		this->m_maskRatio = ratio;
	}
}

/** read the mask image into memory. */
void ReadMaskImage( std::string maskFileName )
{
	typedef itk::ImageFileReader< MaskImageType >					MaskImageReaderType;
	typename MaskImageReaderType::Pointer	maskreader = MaskImageReaderType::New();
	
	maskreader->SetFileName( maskFileName );

	try{
		maskreader->Update();
	}
	catch( itk::ExceptionObject &err ){
		std::stringstream msg("");
		msg<<"An error was caught reading the mask image. Error message:"<<std::endl;
		msg<<err<<std::endl;
		this->WriteToLogFile( msg.str() );
		std::abort();
	}
	
	this->SetMaskImage( maskreader->GetOutput() );
}

/** return the mask image label object. */
LabelMapPointer GetLabelMap()
{
	return this->m_labelMap;
}

/** Inspects the mask image for regions. Takes the largest region and 
 * places it into the label mask for definition of the output image. */
void CalculateMaskImageLabel()
{
	typedef itk::BinaryImageToShapeLabelMapFilter< MaskImageType >	BinaryImageToShapeLabelMapFilterType;
	typename BinaryImageToShapeLabelMapFilterType::Pointer binaryToLabelFilter = BinaryImageToShapeLabelMapFilterType::New();
	binaryToLabelFilter->SetInput( this->GetMaskImage() );
	try{
		binaryToLabelFilter->Update();
	}
	catch( itk::ExceptionObject &err ){
		std::stringstream msg("");
		msg<<"An error was caught when calculating the mask image labels.  Error Message:"<<std::endl;
		msg<<err<<std::endl;
		this->WriteToLogfile( msg.str() );
		std::abort();
	}
	
	int labelSizes[ binaryToLabelFilter->GetOutput()->GetNumberOfLabelObjects() ];
	
	for( int i = 0; i << binaryToLabelFilter->GetOutput()->GetNumberOfLabelObjects(); ++i){
		LabelObjectPointer labelObject = binaryToLabelFilter->GetOutput()->GetNthLabelObject(i);
		labelSizes[i] = labelObject->Size();
	}
	
	int indexLargestLabel = *std::max_element( labelSizes, labelSizes+binaryToLabelFilter->GetOutput()->GetNumberOfLabelObjects() );
	
	LabelMapPointer labelMap = LabelMapType::New();
	labelMap->SetOrigin( this->GetMaskImage()->GetOrigin() );
	labelMap->SetSpacing( this->GetMaskImage()->GetSpacing() );
	labelMap->SetRegions( this->GetMaskImage()->GetLargestPosssibleRegion() );
	labelMap->AddLabelObject( binaryToLabelFilter->GetOutput()->GetNthLabelObject( indexLargestLabel ) );
		
	this->SetLabelMap( labelMap );
}

/** Set the mask image label object. This object can be used to find the
 * center of the labeled region, and its bounding box. Returns true if
 * the labeled region was placed into or was already in m_maskLabel. 
 * Returns false otherwise. */
bool SetLabelMap( LabelMapPointer labelMap )
{
	if( this->m_labelMap.GetPointer() != labelMap ){
		this->m_labelMap = labelMap;
		return true;
	}
	else if( this->m_labelMap.GetPointer() == labelMap){
		return true;
	}
	else{
		return false;
	}
}

/** Set the mask image. Returns true if the provided image was set as or
 * already was stored in m_maskImage. Returns false if the setting was 
 * not performed for any reason. */
bool SetMaskImage( MaskImagePointer maskImage )
{
	if( this->m_maskImage.GetPointer() != maskImage ){
		this->m_maskImage = maskImage;
		return true;
	}
	else if( this->m_maskImage.GetPointer() == maskImage){
		return true;
	}
	else{
		return false;
	}
}

MaskImagePointer GetMaskImage()
{
	return this->m_MaskImage;
}

private:

DataImagePointer			m_DataImage;
double						m_errorRadius;
double						m_displacementErrorTolerance;
double						m_strainErrorTolerance;
double						m_maxMeticValue; // this because I'm using the normalized x-correlation coefficient metric
vtkSmartPointer<vtkIdList>	m_pointsList;
MaskImagePointer			m_maskImage;
RegionOverlapType			m_regionOverlap;
MaskRatioType				m_maskRatio;
LabelMapPointer				m_labelMap;



};

#endif

