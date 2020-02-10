# Refinement criteria

## Introduction
There are a range of hard-coded refinement options which you can choose between. 

As refinement criteria tend to be tailored for different problems, it's quite possible that unfortunately the existing options will not quite do what you want. It shouldn't be too difficult to add your own though - if you need to do this and need some help then please get in touch.

# In built options
Here is a quick summary of what currently exists/where/how to find it.

The important code is found in `AMRLevelMushyLayerRegrid::tagCells()`. 

The refinement method is determined by `(m_)opt.refinementMethod` and read in MushyLayerSubcycleUtils::getAMRFactory(). Unfortunately it's not as simple as something like main.refineMethod=... due to backward compatibility issues. I think the following list of inputs options/refinement methods is correct, but I haven't tested it all so apologies if not. In each case, the code will only ever use one of these methods (not combinations of them).

* `main.vel_refine_thresh=X` : refine wherever fluid speed exceeds `X`
* `regrid.tag_channels=1` : try and refine where there are channels according to some hard coded criteria in AMRLevelMushyLayerRegrid::tagCells()
* `regrid.plume_vel=X` and `regrid.plume_salinity=Y` : refine where velocity exceeds `X` and salinity exceeds `Y`, which can be a good indicator of where channels exists if X and Y are chosen appropriately
* `main.taggingVar=X` and `main.refine_thresh=Y`: refine where the undivided gradient of some scalar field exceeds `Y`, where X is the index of the scalar field (0=enthalpy, 1=bulk concentration, 2=temperature... according to the order they are defined in ScalarVars - see [https://amr-softball.github.io/doc/html/mushy_layer_opt_8h.html#afcada9fb65a998951da882b5c10191fe])
* `main.taggingVectorVar=X` and `main.refine_thresh=Y`: like (4) but for vector fields.
* `regrid.tag_mush_channels=1` : refine wherever porosity < 1 on level 1, then try and refine around channels on higher levels according to some hardcoded criteria.
* `regrid.tag_channels_composite=1` : yet another hard coded attempt to refine around channels/mushy regions. Here, we use a single criteria on all levels. 

Each of these options should be mapped to one of the RefinementMethod options listed here: [https://amr-softball.github.io/doc/html/mushy_layer_opt_8h.html#ad767f132740c90bc3329dc12556cd562] within `MushyLayerSubcycleUtils::getAMRFactory()`. There are also a few other extra options which are used in some (by no means all) of these criteria, such as

``` 
regrid.tagMushLiquidBoundary
regrid.tagDomainBoundary
regrid.onlyTagPorousCells
regrid.porousCellsShrink
```

I think these are mostly self explanatory. The last one defines an integer number of cells by which the set of cells that has been tagged for refinement is shrunk and then grown again to stop refining around small artifacts. E.g. if regrid.porousCellsShrink=3 and two adjacent random cells were tagged somewhere in the domain, after shrinking the set of tagged cells by 3 this group would have disappeared. 
