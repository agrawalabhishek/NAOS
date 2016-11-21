# Copyright (c) <year> <author> (<email>)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

# Set project source files.
set(SRC
  "${SRC_PATH}/ellipsoidSurfacePoints.cpp"
  "${SRC_PATH}/cubicRoot.cpp"
  "${SRC_PATH}/ellipsoidGravitationalAcceleration.cpp"
  "${SRC_PATH}/computeEllipsoidSurfaceGravitationalAcceleration.cpp"
  "${SRC_PATH}/executeOrbiterAroundURE.cpp"
  "${SRC_PATH}/gslIntegratorOrbiterAroundURE.cpp"
  "${SRC_PATH}/gslIntegratorOrbiterAroundUREPointMassGravity.cpp"
  "${SRC_PATH}/regolithTrajectoryCalculator.cpp"
  "${SRC_PATH}/postAnalysis.cpp"
  "${SRC_PATH}/boostIntegratorRestrictedTwoBodyProblem.cpp"
  "${SRC_PATH}/particleAroundSpheroidAndEllipsoidGravitationalPotential.cpp"
)

# Set project main file.
set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

# Set project test source files.
set(TEST_SRC
  "${TEST_SRC_PATH}/testNaos.cpp"
  "${TEST_SRC_PATH}/testConvertKeplerToCartesian.cpp"
  "${TEST_SRC_PATH}/testConvertCartesianToKepler.cpp"
  "${TEST_SRC_PATH}/testCrossProduct.cpp"
  "${TEST_SRC_PATH}/testMaxRealCubicRoot.cpp"
  "${TEST_SRC_PATH}/testDotProduct.cpp"
  "${TEST_SRC_PATH}/testVectorNorm.cpp"
  "${TEST_SRC_PATH}/testMatrixMultiplication.cpp"
  "${TEST_SRC_PATH}/testMatrixTranspose.cpp"
)
