# Run all convergence tests for the methods paper
from testDiffusiveSolidification import testDiffusiveSolidification
from testUniformPorousConvection import testUniformPorousConvection
from testFixedPorousHole import testFixedPorousHole
from testPorousMushyHole import testPorousMushyHole
from testFixedChillHeleShaw import testFixedChillHeleShaw

testDiffusiveSolidification()
testUniformPorousConvection()
testFixedPorousHole()
testPorousMushyHole()
testFixedChillHeleShaw()


# 5) Fixed chill Hele-Shaw cell


