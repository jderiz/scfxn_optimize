
    # if optimizer in ["forest", "all"]:
    #     print("FOREST")
    #     forest_result = forest_minimize(
    #         func=objective,
    #         dimensions=dimensions,
    #         acq_func="gp_hedge",  # expected Improvement (evaluate whats best to use)
    #         n_calls=n_calls,
    #         x0=default_parameters,
    #         n_jobs=-1,
    #         xi=xi,
    #         kappa=kappa,
    #         callback=[timer_callback, forest_check],
    #     )  # All availiable Cores if aquisition =lbfgs
    #     with open(
    #         "results/res_{}_{}_similar.pkl".format("forest", str(n_calls)), "wb"
    #     ) as file:
    #         pickle.dump(forest_result, file)

    # if optimizer in ["gbrt", "all"]:
    #     print("GradientBoostedTrees")
    #     gbrt_res = gbrt_minimize(
    #         func=objective,
    #         dimensions=dimensions,
    #         acq_func="gp_hedge",
    #         n_calls=n_calls,
    #         x0=default_parameters,
    #         n_jobs=-1,
    #         xi=xi,
    #         kappa=kappa,
    #         callback=[timer_callback, gbrt_check],
    #     )
    #     with open(
    #         "results/res_{}_{}_similar.pkl".format("gbrt", str(n_calls)), "wb"
    #     ) as file:
    #         pickle.dump(gbrt_res, file)

    # if optimizer in ["gp", "all"]:
    #     print("Gaussian Processes")
    #     gp_res = gp_minimize(
    #         func=objective,
    #         dimensions=dimensions,
    #         acq_func="gp_hedge",
    #         n_calls=n_calls,
    #         n_jobs=-1,
    #         x0=default_parameters,
    #         xi=xi,
    #         kappa=kappa,
    #         callback=[timer_callback, gp_check],
    #     )
    #     with open(
    #         "results/res_{}_{}_similar.pkl".format("gp", str(n_calls)), "wb"
    #     ) as file:
    #         pickle.dump(gp_res, file)
