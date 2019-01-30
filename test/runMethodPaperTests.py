# Run all convergence tests for the methods paper
from testDiffusiveSolidification import testDiffusiveSolidification
from testUniformPorousConvection import test_uniform_porous_convection
from testFixedPorousHole import testFixedPorousHole
from testPorousMushyHole import testPorousMushyHole
from testHeleShawFixedChill import testHeleShawFixedChill

testDiffusiveSolidification()
test_uniform_porous_convection([])
testFixedPorousHole()
testPorousMushyHole([])
testHeleShawFixedChill([])


# 5) Fixed chill Hele-Shaw cell


