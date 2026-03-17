"""
Potential Flow Analyzer
=======================
A 2D potential flow visualization tool that allows the user to superimpose
classical fluid flow elements to model and visualize complex flow fields.

Supported flow elements:
    S1 - Uniform Flow
    S2 - Source
    S3 - Sink
    S4 - Doublet
    S5 - Vortex

The program computes the combined stream function symbolically, evaluates
it numerically on a 2D grid, detects the body surface, and plots the
resulting streamlines, stagnation points, and flow element origins.

Dependencies:
    numpy, matplotlib, sympy

Author: Aadi Tripathi
Date: March 2025
"""



# ── External Libraries ──────────
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp



"""
    Below builds a symbolic stream function for a given flow element type.

    Args:
        function_name (str): The flow element type. One of 'S1', 'S2', 'S3', 'S4', 'S5'.
        index (int): The position of this function in the summation order.
        U (float, optional): Uniform flow velocity. Required when building S2 on top of S1.

    Returns:
        S (sympy expression): The symbolic stream function.
        u (sympy expression): Velocity component in the x-direction (dS/dy).
        v (sympy expression): Velocity component in the y-direction (-dS/dx).
        dict: Origin coordinates of this flow element {a_index: a, b_index: b}.
        Returns (None, None, None, None) if input is invalid.
    """



def create_s_function(function_name, index, U=None):
   # Define symbolic variables for x and y
   x, y = sp.symbols('x y')


   # Get the coordinates of the origin for the current function
   try:
       a = float(input(f"Enter x-coordinate of origin for {function_name}: "))
       b = float(input(f"Enter y-coordinate of origin for {function_name}: "))
   except ValueError:
       print("Invalid input. Please enter numeric values.")
       return None, None, None, None


    # Shift coordinates so polar conversion is centered on this element's origin
    # instead of the global origin (0, 0)
   x_shifted = x - a
   y_shifted = y - b


   # Calculate the radial distance and angle
   r = sp.sqrt(x_shifted ** 2 + y_shifted ** 2)
   theta = sp.atan2(y_shifted, x_shifted)


   # ── S1: Uniform Flow ──────────────────────────────────────────────────────
    # Models a constant velocity field flowing horizontally across the whole plane.
    # This is almost always the base layer that other elements are added on top of.
    # Formula: ψ = U * r * sin(θ)  —  equivalent to ψ = U * y in Cartesian form
   if function_name == 'S1':
       try:
           U = float(input("Enter uniform flow velocity (U): "))
       except ValueError:
           print("Invalid input. Please enter a numeric value.")
           return None, None, None, None
       # Uniform flow stream function
       S = U * r * sp.sin(theta)
   
   # ── S2: Source ────────────────────────────────────────────────────────────
    # Models fluid radiating outward equally in all directions from a single point.
    # Q is the source strength — how much fluid is being emitted per unit time.
    # When combined with S1 (uniform flow) this produces the Rankine half-body,
    # which models flow around a blunt nosed shape like the front of a submarine.
    # Formula: ψ = (Q * θ) / 2π
   elif function_name == 'S2':
       if U is None:
           print("Error: U must be provided for S2.")
           return None, None, None, None
       try:
           Q = float(input("Enter source strength (Q): "))
       except ValueError:
           print("Invalid input. Please enter numeric values.")
           return None, None, None, None
       # Calculate stagnation distance b — the distance from the source to the
        # stagnation point where the uniform flow and source flow cancel out
       b = Q / (2 * np.pi * U)  
       S = (Q * theta) / (2 * sp.pi)
    
    # ── S3: Sink ──────────────────────────────────────────────────────────────
    # The opposite of a source — fluid converges into a single point from all
    # directions. Mathematically identical to S2 but with a negative sign.
    # Q is the sink strength — how much fluid is being absorbed per unit time.
    # A source and sink placed close together begin to approximate a doublet (S4).
    # Formula: ψ = -(Q * θ) / 2π
   elif function_name == 'S3':
       try:
           Q = float(input(f"Enter value for Q for {function_name}: "))
       except ValueError:
           print("Invalid input. Please enter a numeric value.")
           return None, None, None, None
        # Sink stream function — negative of source
       S = (-Q / (2 * sp.pi)) * theta
    
    # ── S4: Doublet ───────────────────────────────────────────────────────────
    # Models the limiting case of a source and sink pushed infinitely close together
    # with their strengths increasing to compensate. The result is a purely radial
    # flow pattern that decays with distance from the origin.
    # k is the doublet strength.
    # Key use case: S1 + S4 combined produces flow around a perfect cylinder.
    # The sin(θ)/r term is what causes the strength to decay with distance —
    # strong near the center, negligible far away.
    # Formula: ψ = -(k * sin(θ)) / (2π * r)
   elif function_name == 'S4':
       try:
           k = float(input(f"Enter value for k for {function_name}: "))
       except ValueError:
           print("Invalid input. Please enter a numeric value.")
           return None, None, None, None
       # Doublet stream function
       S = (-k / (2 * np.pi * r)) * sp.sin(theta)
   
   # ── S5: Vortex ────────────────────────────────────────────────────────────
    # Models fluid rotating around a central point in concentric circles.
    # R is the reference radius and w is the angular rotation speed.
    # The log(r/R) term means the vortex effect spreads across the whole plane
    # but weakens with distance — unlike a doublet it never fully disappears.
    # Key use case: S1 + S4 + S5 combined produces a rotating cylinder which
    # generates lift — this is the fundamental principle behind how airfoils work
    # and is known as the Kutta-Joukowski theorem.
    # Formula: ψ = (-2 * R² * w) * log(r/R)
   elif function_name == 'S5':
       try:
           R = float(input(f"Enter value for R for {function_name}: "))
           w = float(input(f"Enter value for w for {function_name}: "))
       except ValueError:
           print("Invalid input. Please enter numeric values.")
           return None, None, None, None
       # Vortex stream function
       S = ((-4 * sp.pi * R ** 2 * w) / (2 * sp.pi)) * sp.log(r / R)
   else:
       print(f"Invalid function name: {function_name}")
       return None, None, None, None


   # ── Velocity Components ───────────────────────────────────────────────────
    # In potential flow theory velocity components are derived directly from the
    # stream function by differentiation. This is done symbolically here so that
    # the exact analytical derivatives are used rather than numerical approximations.
    # u =  dψ/dy  — velocity in the x-direction
    # v = -dψ/dx  — velocity in the y-direction
    # The negative sign on v comes from the definition of the stream function.
   u = sp.diff(S, y)  # Velocity component in the y-direction
   v = -sp.diff(S, x)  # Velocity component in the x-direction


   # Return the symbolic stream function, both velocity components, and the
    # origin coordinates of this element so they can be plotted later in main()
   return S, u, v, {f'a_{index}': a, f'b_{index}': b}




# Function to detect wall streamlines (zero contour of the stream function)
"""
    Detects the body surface (wall) from the computed stream function grid.

    The wall is defined as the zero contour of the stream function (ψ=0).
    Points with near-zero gradients are also flagged as they represent
    stagnant interior regions where plotting streamlines is not meaningful.

    Args:
        psi_values (ndarray): 2D grid of stream function values (300x300).

    Returns:
        wall_mask (ndarray): Boolean 2D grid where True indicates a wall
                             or interior point that should not be plotted.
    """
def detect_wall_streamline(psi_values):
   
   # Find all points where ψ is close to zero — this is the body surface.
    # atol=1e-3 is a small tolerance since floating point math never gives
    # exactly 0.0, so anything within 0.001 is treated as the wall.
   zero_contour = np.isclose(psi_values, 0, atol=1e-3)  
   
    # Compute how fast ψ is changing across the grid in both directions.
    # Think of this as finding the slope of the stream function at every point.
   grad_x, grad_y = np.gradient(psi_values)  
   
    # Combine both gradient directions into a single magnitude value per point
   grad_magnitude = np.sqrt(grad_x ** 2 + grad_y ** 2)  
   
   # Flag points where ψ is barely changing — these are stagnant interior
    # regions inside the body where streamlines should not be drawn
   low_gradient = grad_magnitude < 1e-3  
   
   # Final wall mask combines both conditions — a point is considered wall
    # if it is ON the ψ=0 contour OR inside a stagnant interior region
   wall_mask = zero_contour | low_gradient  # Mark wall streamlines
   
   return wall_mask




# Function to compute psi values at given X, Y mesh points
"""
    Evaluates the total stream function numerically across the entire 2D grid.

    This is the moment where the symbolic math transitions into actual numbers.
    NumPy evaluates all 90,000 grid points simultaneously making it very fast.

    Args:
        X (ndarray): 2D grid of x coordinates (300x300) from np.meshgrid.
        Y (ndarray): 2D grid of y coordinates (300x300) from np.meshgrid.
        sum_function_np (callable): The fully assembled numerical stream function
                                    converted from symbolic form by safe_lambdify.

    Returns:
        psi_values (ndarray): 2D grid of stream function values at every point.
    """
def compute_psi_values(X, Y, sum_function_np):
   # Evaluate the stream function at every point on the grid simultaneously.
    # This single line replaces what would otherwise be a nested loop over
    # all 300x300 = 90,000 grid points — NumPy handles it in one vectorized call.
   return sum_function_np(X, Y)




# Function to plot streamlines outside the wall (masked regions are NaN)
"""
    Plots streamlines across the 2D grid while masking the body interior.

    Uses NaN replacement to prevent matplotlib from drawing streamlines
    inside the body surface. Matplotlib's streamplot automatically skips
    any grid point that contains a NaN value.

    Args:
        ax (matplotlib Axes): The figure panel to draw streamlines onto.
        X (ndarray): 2D grid of x coordinates (300x300).
        Y (ndarray): 2D grid of y coordinates (300x300).
        U (ndarray): 2D grid of x-direction velocity values (300x300).
        V (ndarray): 2D grid of y-direction velocity values (300x300).
        wall_mask (ndarray): Boolean grid where True marks interior/wall points.

    Returns:
        None. Draws directly onto the provided matplotlib Axes object.
    """
def plot_streamlines_outside_wall(ax, X, Y, U, V, wall_mask):
   # Replace velocity values with NaN at every wall or interior point.
    # np.where reads: "if wall_mask is True at this point use NaN,
    # otherwise keep the real velocity value."
    # matplotlib's streamplot silently skips NaN points so nothing gets
    # drawn inside the body without needing any extra logic.
   U_masked = np.where(wall_mask, np.nan, U)  # Mask U values at the wall
   V_masked = np.where(wall_mask, np.nan, V)  # Mask V values at the wall
   
   # Draw streamlines across the masked velocity field.
    # density=2.5 controls how many lines are drawn — higher means more
    # lines and more detail but can look cluttered above 3.0
   ax.streamplot(X, Y, U_masked, V_masked, color='black', density=2.5)  # Plot streamlines




# Safely lambdify the symbolic expression for numerical evaluation
"""
    Safely converts a SymPy symbolic expression into a fast NumPy numerical function.

    This is the bridge between the symbolic algebra world (SymPy) and the
    numerical computation world (NumPy). Without this conversion the stream
    function cannot be evaluated across the grid.

    Think of it as translating from algebra language into calculator language.

    Args:
        symbols (tuple): The symbolic variables to convert, typically (x, y).
        expression (sympy expression): The symbolic expression to convert.
        modules (str): The numerical backend to use. Defaults to 'numpy'.

    Returns:
        callable: A fast NumPy function that can be evaluated across the grid.
                  If conversion fails returns a zero function to prevent crashing.
    """
def safe_lambdify(symbols, expression, modules='numpy'):
   try:
       # Convert the symbolic expression to a numerical function.
        # sp.lambdify is SymPy's built in tool for this conversion —
        # the result behaves like a normal Python function but runs at NumPy speed.
       return sp.lambdify(symbols, expression, modules)
   except Exception as e:
       print(f"Lambdify error: {e}")
       # If conversion fails for any reason return a function that gives zeros
        # everywhere — this prevents the whole program from crashing and allows
        # the rest of the visualization to still run with a blank field
       return lambda x, y: np.zeros_like(x)  # Return a zero function if lambdify fails




# Main function that drives the flow analysis
"""
    Main driver function for the Potential Flow Analyzer.

    Handles all user input, builds the total stream function by superimposing
    the selected flow elements, computes the flow field numerically, detects
    the body surface, and generates the final visualization.
    """
def main():
    # ── Initialization ────────────────────────────────────────────────────────
    # All accumulators start at zero — each selected flow element gets added
    # on top symbolically, building up the total stream function one piece at a time
   sum_function = 0          # Accumulates the total stream function ψ
   sum_u = 0                 # Accumulates the total x-direction velocity
   sum_v = 0                 # Accumulates the total y-direction velocity
   selected_functions = []   # Tracks which S functions were chosen, used for plot title
   variables = {}            # Stores origin coordinates of each element for plotting


   # Ask how many functions the user wants to sum
   num_functions = int(input("Enter the number of S functions to sum: "))


   rankine_params = {}  # Dictionary to store parameters for Rankine flow

 # ── Input Loop ────────────────────────────────────────────────────────────
    # Loop through each function the user wants to add.
    # The while True loop keeps asking until a valid S function name is entered.
   for i in range(num_functions):
       while True:
           
           # Ask the user to input which stream function they want to add
           s = input(f"Enter S{i + 1} (S1, S2, S3, S4, or S5): ").upper()
           if s in ['S1', 'S2', 'S3', 'S4', 'S5']:
               
               # Store U and Q in rankine_params so they are available for
                # the Rankine half-body analysis printed later in main()
               if s == 'S1':
                   rankine_params['U'] = float(input("Enter uniform flow velocity (U): "))
               elif s == 'S2':
                   # S2 requires S1 to be defined first since it needs U
                   if 'U' not in rankine_params:
                       print("Error: Please define S1 with U before using S2.")
                       continue
                   rankine_params['Q'] = float(input("Enter source strength (Q): "))
                   # Stagnation distance — distance from source to stagnation point
                   rankine_params['b'] = rankine_params['Q'] / (2 * np.pi * rankine_params['U'])


               # Build this element's stream function and velocity components
                # then accumulate them into the running symbolic total
               S, u, v, func_vars = create_s_function(s, i, U=rankine_params.get('U'))
               if S is None:
                   continue


               sum_function += S  # Add this element's stream function to total
               sum_u += u  # Add this element's x-velocity to total
               sum_v += v   # Add this element's y-velocity to total
               variables.update(func_vars)      # Store origin coordinates
               selected_functions.append(s)      # Record which function was added
               break
           else:
               print("Invalid input. Please enter S1, S2, S3, S4, or S5.")


    # ── Rankine Half-Body Analysis ────────────────────────────────────────────
    # If the user selected exactly S1 + S2, print the analytical Rankine results.
    # This is a special case because S1 + S2 has a known closed form solution
    # for the stagnation point and body outline equation.
   if set(selected_functions) == {'S1', 'S2'}:
       print("\n--- Rankine Half-Body Flow Analysis ---")
       U, Q, b = rankine_params['U'], rankine_params['Q'], rankine_params['b']


       # Stagnation velocity — the uniform flow speed at the stagnation point.
        # At this point the source flow exactly cancels the uniform flow giving
        # zero net velocity. This is where the body surface begins.
       U_stagnation = Q / (2 * np.pi * b)
       print(f"Stagnation Velocity: U_stagnation = {U_stagnation}")
       print(f"Stagnation Point: (x, y) at r=-{b}, θ=π")


       # The body outline equation in polar form — describes the shape of the
        # Rankine half-body surface at every angle θ
       print("\nBody Outline Equation:")
       print("r = b(π - θ) / sin(θ)")
       print(f"Where b={b}")


   # Print the full symbolic stream function to the terminal so the user
    # can verify the superposition looks correct before the plot appears
   print("\nFinal function as the sum of selected S's:")
   print(sp.pretty(sum_function))


   # ── Symbolic to Numerical Conversion ─────────────────────────────────────
    # Convert all three symbolic expressions to fast NumPy numerical functions.
    # This must be done before any grid evaluation can happen
   x, y = sp.symbols('x y')
   sum_function_np = safe_lambdify((x, y), sum_function)  # Stream function ψ
   u_np = safe_lambdify((x, y), sum_u)           # x-direction velocity
   v_np = safe_lambdify((x, y), sum_v)             # y-direction velocity


   # ── Grid Setup ────────────────────────────────────────────────────────────
    # Build a 300x300 grid covering -5 to 5 in both x and y directions.
    # np.meshgrid converts two 1D arrays into two 2D coordinate grids so that
    # the stream function can be evaluated at every point simultaneously
   x_vals = np.linspace(-5, 5, 300)
   y_vals = np.linspace(-5, 5, 300)
   X, Y = np.meshgrid(x_vals, y_vals)


   # ── Grid Evaluation ───────────────────────────────────────────────────────
    # Evaluate the stream function and velocity components across the full grid.
    # Wrapped in try/except so that if evaluation fails the program still runs
    # and produces a blank plot rather than crashing entirely.
   try:
       psi_values = compute_psi_values(X, Y, sum_function_np)  # ψ at every point
       U = u_np(X, Y)  # x-velocity at every point
       V = v_np(X, Y)  # y-velocity at every point
   except Exception as e:
       print(f"Error computing values: {e}")   # Fall back to zero arrays so the rest of the program can still run
       psi_values = np.zeros_like(X)
       U = np.zeros_like(X)
       V = np.zeros_like(X)


   # Detect the body surface from the computed ψ grid
   wall_mask = detect_wall_streamline(psi_values)


   # ── Figure Setup ──────────────────────────────────────────────────────────
    # Create a 12x12 inch square figure — equal dimensions keep the aspect ratio
    # correct so circular flow patterns don't appear stretched or squished
   fig, ax = plt.subplots(figsize=(12, 12))


    # ── Streamlines ───────────────────────────────────────────────────────────
    # Draw streamlines across the velocity field, masking the body interior.
    # Wrapped in try/except so a plotting failure doesn't crash the whole program.
   try:
       # Plot streamlines outside the wall
       plot_streamlines_outside_wall(ax, X, Y, U, V, wall_mask)
   except Exception as e:
       print(f"Streamline plotting error: {e}")


    # ── Body Surface ──────────────────────────────────────────────────────────
    # Draw the body outline as a red contour line at ψ=0.
    # The [0] argument tells matplotlib to draw only the zero contour —
    # this is the streamline that defines the body surface shape.
   ax.contour(X, Y, psi_values, [0], colors='red', linewidths=2)


   # ── Stagnation Points ─────────────────────────────────────────────────────
    # Find stagnation points by looking for where BOTH velocity components are
    # near zero simultaneously — meaning the flow has come to a complete stop.
    # Only points outside the wall are valid stagnation points.
   stagnation_mask = (np.abs(U) < 1e-3) & (np.abs(V) < 1e-3)
   valid_stagnation_mask = stagnation_mask & ~wall_mask
   stagnation_points = np.column_stack((X[valid_stagnation_mask], Y[valid_stagnation_mask]))


   # Plot stagnation points as blue dots with labels
   if len(stagnation_points) > 0:
       ax.scatter(stagnation_points[:, 0], stagnation_points[:, 1], color='blue', s=50,
                  zorder=5, label='Stagnation Points')


       # Annotate each stagnation point with its number
        # zorder=5 ensures the dots appear on top of the streamlines
       for i, (x_val, y_val) in enumerate(stagnation_points):
           ax.annotate(f'Stagnation Point {i + 1}', (x_val, y_val),
                       xytext=(5, 5), textcoords='offset points', color='blue')
   else:
       print("No stagnation points found.")


   # ── Flow Element Origins ──────────────────────────────────────────────────
    # Plot the origin point of each flow element as a green dot so the user
    # can see exactly where each S function was placed in the flow field
   try:
       for i, s in enumerate(selected_functions):
           a = variables[f'a_{i}']
           b = variables[f'b_{i}']
           ax.scatter(a, b, color='green', s=100, zorder=4)
           ax.annotate(f"{s} ({a},{b})", (a, b),
                       xytext=(5, -10), textcoords='offset points')
   except Exception as e:
       print(f"Error marking function origins: {e}")


   # ── Plot Formatting ───────────────────────────────────────────────────────
    # Standard axis labels, title built from selected function names,
    # equal aspect ratio to preserve circular flow geometry,
    # and a light grid for readability
   ax.set_xlabel('X', fontsize=12)
   ax.set_ylabel('Y', fontsize=12)
   ax.set_title(f"Streamlines and Stagnation Points: {' + '.join(selected_functions)}", fontsize=14)
   ax.set_xlim(-5, 5)
   ax.set_ylim(-5, 5)
   ax.set_aspect('equal')
   ax.grid(True, linestyle='--', alpha=0.7)


    # Render the final plot
   plt.tight_layout()
   plt.show()




# ── Entry Point ───────────────────────────────────────────────────────────────
# Standard Python entry point check — ensures main() only runs when this file
# is executed directly, not when it is imported as a module by another script
if __name__ == "__main__":
   main()
