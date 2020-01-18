# Run all convergence tests for the methods paper
from testDiffusiveSolidification import test_diffusive_solidification
from testUniformPorousConvection import test_uniform_porous_convection
from testFixedPorousHole import test_fixed_porous_hole
from testPorousMushyHole import test_porous_mushy_hole
from testHeleShawFixedChill import test_hele_shaw_fixed_chill

test_diffusive_solidification()
test_uniform_porous_convection([])
test_fixed_porous_hole()
test_porous_mushy_hole([])
test_hele_shaw_fixed_chill([])


# 5) Fixed chill Hele-Shaw cell


