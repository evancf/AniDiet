# AniDiet

The broad goal of this project is to quantify global animal food webs, aiming to understand macroecological variation in food web properties and how global change drivers affect food webs. Our approach is to pair data on animal diets assembled from the scientific literature with machine learning models to build food webs at any site around the world, thereby linking species composition to ecosystem functioning.

The first aim was to reconstruct terrestrial mammal food webs since the last interglacial, and we have developed a manuscript "Collapse of terrestrial mammal food webs since the Late Pleistocene" based on the resulting analyses. This repository includes numbered R scripts. Scripts 01-07 provide the reproducible workflow to develop the datasets analyzed in the paper (scripts 04 and 05 show another direction to generate food webs based on observed interactions only, without any trait-based modeling of predator-prey interactions). Script 08 shows the deep learning model to generate food webs given any mammal species assemblage and uses this model to hindcast how food webs would be structured had they been unaffected by extinction since the Late Pleistocene. Script 09 calculates metrics to quantify how food webs have changed over time due to extinction and range changes, and how they would be affected by endangered species extinction.


