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
#include <vtkStructuredGrid.h>
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

template<typename TFixedImage, typename TMovingImage, typename TMaskImage>
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
typedef typename DIC<TFixedImage,TMovingImage>::TransformInitializerType	TransformInitializerType;
typedef typename DIC<TFixedImage,TMovingImage>::TransformType				TransformType;

typedef			 TMaskImage							MaskImageType;
typedef	typename MaskImageType::Pointer				MaskImagePointer;
typedef	typename MaskImageType::ConstPointer		MaskImageConstPointer;

typedef		vtkSmartPointer<vtkStructuredGrid>		DataImagePointer;
typedef		vtkSmartPointer<vtkPoints>				DataImagePointsPointer;
typedef		vtkSmartPointer<vtkDoubleArray>			DataImagePixelPointer;

typedef		double									RegionOverlapType;
typedef		double									MaskRatioType;


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
	m_maskImage = 0; // Must be provided by the user
	m_regionOverlap = 0;
	m_maskRatio = 0.3; // the minimum value of region/mask overlap to evaluate the region
}

void CreateEmptyDataImage()
{
	// Find the middle of the image.
	// Find the number of Regions in each direction in the image
	//		use m_regionOverlap
	// create vtkStructuredGrid to with the correct bounds and origin
}

void CalculateInitialMovingImageRgionList()
{
	// loop through the points in the data image
	// check if the region overlaps the mask image
	//		use m_maskRatio
	// if it does not overlap the mask use vtkStructuredGrid::BlankPoint to turn off the point
	// it is does overlap, add the region, defined by m_IRRadius to the region list
}

void CalculateInitialFixedImageRegionList()
{
	// loop through the points in the data image
	// check if the region overlaps the mask image
	//		use m_maskRatio
	// if it does not overlap the mask use vtkStructuredGrid::BlankPoint to turn off the point
	// it is does overlap, add the region, defined by m_IRRadius to the region list
}

FixedImageRegionType GetGlobalRegistrationRegion()
{
	// return the bounding region of the m_DataImage
}

void SetMaskRatio(double ratio)
{
	if( ration != m_maskRatio ){
		this->m_maskRatio = ratio;
	}
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


}

#endif

