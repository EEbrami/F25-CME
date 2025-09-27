# Example Monte Carlo Simulation for Computational Economics
# This script demonstrates Monte Carlo methods for economic analysis
# Course: Econ-81360, Fall 2025

using Random, Distributions, Statistics, Plots, LinearAlgebra

"""
Monte Carlo simulation of a simple economic model
"""

# Set random seed for reproducibility
Random.seed!(42)

function monte_carlo_portfolio(n_simulations=10000, n_periods=252)
    """
    Monte Carlo simulation of portfolio returns
    
    Parameters:
    - n_simulations: Number of simulation runs
    - n_periods: Number of time periods (e.g., trading days in a year)
    
    Returns:
    - Array of final portfolio values
    """
    
    # Portfolio parameters
    initial_value = 100000  # Initial portfolio value
    expected_return = 0.08  # Annual expected return
    volatility = 0.20       # Annual volatility
    
    # Convert to daily parameters
    daily_return = expected_return / n_periods
    daily_volatility = volatility / sqrt(n_periods)
    
    # Storage for results
    final_values = zeros(n_simulations)
    
    println("Running Monte Carlo simulation...")
    println("Number of simulations: $n_simulations")
    println("Time periods: $n_periods")
    println("Expected annual return: $(expected_return*100)%")
    println("Annual volatility: $(volatility*100)%")
    
    for sim in 1:n_simulations
        portfolio_value = initial_value
        
        for period in 1:n_periods
            # Generate random return
            random_return = rand(Normal(daily_return, daily_volatility))
            
            # Update portfolio value
            portfolio_value *= (1 + random_return)
        end
        
        final_values[sim] = portfolio_value
    end
    
    return final_values
end

function analyze_simulation_results(final_values, initial_value=100000)
    """
    Analyze Monte Carlo simulation results
    """
    
    println("\n" * "="^50)
    println("SIMULATION RESULTS ANALYSIS")
    println("="^50)
    
    # Calculate statistics
    mean_value = mean(final_values)
    median_value = median(final_values)
    std_value = std(final_values)
    min_value = minimum(final_values)
    max_value = maximum(final_values)
    
    # Calculate percentiles
    percentiles = [5, 25, 75, 95]
    percentile_values = [quantile(final_values, p/100) for p in percentiles]
    
    # Print results
    println("Initial Portfolio Value: \$$(round(Int, initial_value))")
    println("\nFinal Portfolio Value Statistics:")
    println("  Mean:     \$$(round(Int, mean_value))")
    println("  Median:   \$$(round(Int, median_value))")
    println("  Std Dev:  \$$(round(Int, std_value))")
    println("  Minimum:  \$$(round(Int, min_value))")
    println("  Maximum:  \$$(round(Int, max_value))")
    
    println("\nPercentiles:")
    for (p, val) in zip(percentiles, percentile_values)
        println("  $(p)th:     \$$(round(Int, val))")
    end
    
    # Risk metrics
    prob_loss = mean(final_values .< initial_value)
    prob_double = mean(final_values .> 2 * initial_value)
    
    println("\nRisk Metrics:")
    println("  Probability of Loss: $(round(prob_loss*100, digits=2))%")
    println("  Probability of Doubling: $(round(prob_double*100, digits=2))%")
    
    return mean_value, std_value, percentile_values
end

function create_simulation_plots(final_values, initial_value=100000)
    """
    Create visualizations of simulation results
    """
    
    # Create plots
    p1 = histogram(final_values, bins=50, alpha=0.7, 
                   title="Distribution of Final Portfolio Values",
                   xlabel="Final Value (\$)", ylabel="Frequency",
                   label="Simulation Results")
    
    # Add vertical lines for statistics
    vline!([initial_value], color=:red, linewidth=2, label="Initial Value")
    vline!([mean(final_values)], color=:green, linewidth=2, label="Mean")
    vline!([median(final_values)], color=:blue, linewidth=2, label="Median")
    
    # Create cumulative probability plot
    sorted_values = sort(final_values)
    n = length(sorted_values)
    probabilities = (1:n) / n
    
    p2 = plot(sorted_values, probabilities, 
              title="Cumulative Probability Distribution",
              xlabel="Portfolio Value (\$)", ylabel="Cumulative Probability",
              linewidth=2, label="CDF")
    
    # Add reference lines
    hline!([0.5], color=:red, linestyle=:dash, label="50th Percentile")
    vline!([initial_value], color=:red, linestyle=:dash, label="Initial Value")
    
    # Combine plots
    combined_plot = plot(p1, p2, layout=(2,1), size=(800, 600))
    
    # Save plot
    savefig(combined_plot, "monte_carlo_simulation.png")
    println("\nPlot saved as 'monte_carlo_simulation.png'")
    
    return combined_plot
end

function economic_policy_simulation()
    """
    Example of using Monte Carlo for economic policy analysis
    """
    
    println("\n" * "="^50)
    println("ECONOMIC POLICY SIMULATION EXAMPLE")
    println("="^50)
    
    # Parameters for policy intervention
    n_simulations = 5000
    baseline_growth = 0.03  # 3% baseline GDP growth
    policy_effect = Normal(0.01, 0.005)  # Policy adds ~1% growth with uncertainty
    
    baseline_outcomes = []
    policy_outcomes = []
    
    for _ in 1:n_simulations
        # Baseline scenario
        baseline_shock = rand(Normal(0, 0.02))
        baseline_gdp = baseline_growth + baseline_shock
        push!(baseline_outcomes, baseline_gdp)
        
        # Policy scenario
        policy_boost = rand(policy_effect)
        policy_gdp = baseline_growth + baseline_shock + policy_boost
        push!(policy_outcomes, policy_gdp)
    end
    
    println("Policy Analysis Results:")
    println("  Baseline mean growth: $(round(mean(baseline_outcomes)*100, digits=2))%")
    println("  Policy mean growth:   $(round(mean(policy_outcomes)*100, digits=2))%")
    println("  Average policy effect: $(round((mean(policy_outcomes) - mean(baseline_outcomes))*100, digits=2))%")
    
    # Probability that policy improves outcomes
    prob_improvement = mean(policy_outcomes .> baseline_outcomes)
    println("  Probability policy improves growth: $(round(prob_improvement*100, digits=1))%")
    
    return baseline_outcomes, policy_outcomes
end

# Main execution
function main()
    println("Monte Carlo Simulation for Economics")
    println("Course: Econ-81360, Fall 2025")
    
    # 1. Portfolio simulation
    final_values = monte_carlo_portfolio(10000, 252)
    
    # 2. Analyze results
    mean_val, std_val, percentiles = analyze_simulation_results(final_values)
    
    # 3. Create visualizations
    plots = create_simulation_plots(final_values)
    
    # 4. Economic policy example
    baseline, policy = economic_policy_simulation()
    
    println("\nSimulation complete!")
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end