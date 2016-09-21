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
)

# Set project main file.
set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

# Set project test source files.
set(TEST_SRC
  "${TEST_SRC_PATH}/testNaos.cpp"
  "${TEST_SRC_PATH}/testConvertKeplerToCartesian.cpp"
  "${TEST_SRC_PATH}/testCrossProduct.cpp"
  "${TEST_SRC_PATH}/testMaxRealCubicRoot.cpp"
  "${TEST_SRC_PATH}/testDotProduct.cpp"
)
