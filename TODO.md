# Thema Masterarbeit
    + "Developing a Framework for optimizing Rosetta Energy Function weights using Bayesian Optimization"
    + Energy Function
    + ML Model structure stable yes/no
    + Can there be any gradient based method ? levinstein
    + PSO / ANTS for Monte-Carlo Guidance


# TODO: 
    - correlationsmatrix für alle ex_ex läufe
    - ex_ex läufe vergleichen
    - Wie soll die Kandidaten Config gewählt werden?
    - Relax Script!!!
    - Gedanken machen wie das ganze gut generalisiert werden kann, sodass ein anwender minimalen config aufwand hat.
     
## Code
    - add time component to loss
    - add cooldown of xi and kappa
    - add SEED für replizierbarkeit
## Analysis
    - show Convergence of Optimization on each loss.
    - Determine 'best' configuration
    - je mehr gute configs ich nehme desto unwahrscheinlicher das sie alle durch zufall entstanden sind. Wenn dann die Configurationen noch näher beieinander legen als zur outroup dann ist das ein weiteres indiz für "echt" gute configs
    - Entfernungsmß für Configs
    - normalisierte distanz anhand der range jedes gewichts.

## Ideas
    - PSO for guiding Monte Carlo
    - Quantum Algo for light vs conductor based.

## Problems
Sinnvoller Loss bei Design
Vergleich von Resultaten
Wenn ich mit prior 1fach teste, und zwar ausschließlich die subgruppe von der ich mir erfolg verspreche verletzt das nicht das Prinzip?
Warum dauern einige runs 60 min mit der gleichen länge und andere nue 30
## papers to read


## Notes
    - Use Ref15 as default params.
    - The Rotamer Packing only linnear scale up to 4 Threads.

# Clustering of results
normalize weight ranges then cluster results by configs and see if there is any significant 
improvement off result metrices in a particular cluster. 
then from the other side aggregate by metrices and check if there is something significant in the best perfoming configurations.

