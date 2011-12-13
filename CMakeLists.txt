cmake_minimum_required(VERSION 2.6)

project( AnalyzeImages)

FIND_PACKAGE(ITK)
FIND_PACKAGE(VTK)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
          "VTK not found. Please set VTK_DIR.")
ENDIF(VTK_FOUND)

ADD_LIBRARY( DIC DIC.cxx )
ADD_LIBRARY( DICMesh DICMesh.cxx )
ADD_LIBRARY( DICGrid DICGrid.cxx )
ADD_EXECUTABLE( AnalyzeImages AnalyzeImages.cxx)
ADD_EXECUTABLE( TestAlgorithm TestAlgorithm.cxx)


TARGET_LINK_LIBRARIES( AnalyzeImages DIC DICMesh DICGrid ITKCommon ITKIO ITKAlgorithms ITKBasicFilters ITKNumerics ITKSpatialObject vtkHybrid )
TARGET_LINK_LIBRARIES( TestAlgorithm DIC DICMesh DICGrid ITKCommon ITKIO ITKAlgorithms ITKBasicFilters ITKNumerics ITKSpatialObject vtkHybrid )

