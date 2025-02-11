from pysat.solvers import Glucose3, Minisat22  # Assuming PySAT's Glucose3 solver is used
from bidict import bidict  # Import the bidict library
from dataclasses import dataclass
import networkx as nx
import time

@dataclass(frozen=True)
class Variable:
    var_type: str  # Type of the variable ('state', 'op', or 'monitor')
    qubit: int = None  # Qubit ID (for state or op variables)
    axis: str = None  # Specifies the axis ('Z' or 'X') for state variables
    col: int = None  # Column (for state or monitor variables)
    timeframe: int = None  # Timeframe (for all variables)
    op_type: str = None  # Operation type (for op variables)
    monitor_type: str = None  # Monitor type (for monitor variables)
    aux_type: str = None  # Auxiliary variable type (for auxiliary variables)
    N: int = None # Number of Qubit to ensure which row the variable it is for state variables

    def format_flag(self) -> str:
        """
        Formats the variable into a human-readable flag for debugging.
        """
        if self.var_type == "state":
            if self.axis == "Z":
                return f"M^{self.timeframe}_{{{self.qubit},{self.col}}} (Z)"
            elif self.axis == "X":
                return f"M^{self.timeframe}_{{{self.qubit + self.N},{self.col}}} (X)"
            else:
                return f"M^{self.timeframe}_{{{self.qubit},{self.col}}}"
        elif self.var_type == "aux":
            return f"{self.aux_type}^{self.timeframe}_{self.col}"
        elif self.var_type == "monitor":
            if self.monitor_type == "completed":
                return f"completed^{self.timeframe}_{self.col}"
            elif self.monitor_type == "property":
                return f"property^{self.timeframe}"
            elif self.monitor_type == "valid":
                return f"valid^{self.timeframe}_{self.col}"
        elif self.var_type == "op":
            if self.op_type in {"S", "H"}:
                return f"{self.op_type}^{self.timeframe}_{self.qubit}"
            elif self.op_type == "CX":
                return f"CX^{self.timeframe}_{{{self.qubit[0]},{self.qubit[1]}}}"
        return "unknown"

    def __str__(self):
        """
        Provides a user-friendly string representation for the variable.
        """
        return self.format_flag()

    def __repr__(self):
        """
        Provides a concise representation for debugging purposes.
        Includes all attributes with non-None values.
        """
        attrs = {key: value for key, value in self.dict__.items() if value is not None}
        attrs_str = ", ".join(f"{key}={value}" for key, value in attrs.items())
        return f"Variable({attrs_str})"

class tableau_resynthesis:
    def __init__(self, tableau):
        # Initialize the SAT formula and solver
        self.cnf = []  # CNF formula as a list of clauses (each clause is a list of literals)
        self.assumption = [] 
        self.solver =  Minisat22()  # SAT solver instance
        self.qubitNum = len(tableau) // 2  # Number of qubits
        self.rotationNum = len(tableau[0]) if len(tableau) > 0 else 0  # Number of rotations
        self.tableau = tableau  # Store the tableau
        self.current_state = None  # Cache for the current state matrix
        self.MaxDepth = 1000  # Maximum depth for unbounded model checking
        self.DG = nx.DiGraph()

        # Use bidirectional mapping for variables
        self.var_bimap = bidict()  # Maps variable indices to Variable objects and vice versa
        self.next_var_index = 1  # Keeps track of the next available SAT variable index

    def build_rotation_dependency_graph(self):
        """
        Builds a dependency graph for a set of Pauli rotations based on their commutation relationships.
        The graph is stored in `self.DG`.

        Returns:
            nx.DiGraph: A directed acyclic graph (DAG) representing dependencies between rotations.
        """
        # Create a directed graph
        graph = nx.DiGraph()

        # Add nodes for each rotation
        graph.add_nodes_from(range(self.rotationNum))

        # Determine dependencies
        for i in range(self.rotationNum):
            for j in range(i + 1, self.rotationNum):
                # Extract Z and X parts of the columns
                z_i = [self.tableau[q][i] for q in range(self.qubitNum)]
                x_i = [self.tableau[q + self.qubitNum][i] for q in range(self.qubitNum)]
                z_j = [self.tableau[q][j] for q in range(self.qubitNum)]
                x_j = [self.tableau[q + self.qubitNum][j] for q in range(self.qubitNum)]

                # Compute symplectic inner product
                symplectic_inner_product = sum((z_i[q] * x_j[q] + x_i[q] * z_j[q]) % 2 for q in range(self.qubitNum)) % 2

                # Add edge if anti-commute
                if symplectic_inner_product == 1:
                    graph.add_edge(i, j)

        # Reduce transitive edges
        self.DG = nx.transitive_reduction(graph)
        return self.DG    
    
    def print_graph(self):
        """
        Prints the dependency graph (self.DG) in a readable format.
        """
        if self.DG is None:
            print("The dependency graph has not been built yet.")
            return

        print("Dependency Graph:")
        print("Nodes (Rotations):")
        for node in self.DG.nodes:
            print(f"  Rotation {node}")

        print("\nEdges (Dependencies):")
        for edge in self.DG.edges:
            print(f"  Rotation {edge[0]} -> Rotation {edge[1]}")

    def print_all_topological_orders(self):
        """
        Prints all possible topological orders of the dependency graph (self.DG).
        """
        if self.DG is None:
            print("The dependency graph has not been built yet.")
            return
        
        # Check if the graph is acyclic
        if not nx.is_directed_acyclic_graph(self.DG):
            print("The graph is not a DAG, so topological sorting is not possible.")
            return
        
        def backtrack(order, visited, stack):
            """
            Recursive helper function to generate all topological orders.
            """
            if len(order) == len(self.DG.nodes):
                # Print the current valid order
                print(order)
                return
            
            # Iterate over nodes with no incoming edges
            for node in list(stack):
                # Choose the current node
                stack.remove(node)
                order.append(node)
                
                # Temporarily remove edges from the graph
                removed_edges = []
                for successor in list(self.DG.successors(node)):
                    removed_edges.append((node, successor))
                    self.DG.remove_edge(node, successor)
                    if self.DG.in_degree(successor) == 0:
                        stack.add(successor)
                
                # Recur with the updated state
                backtrack(order, visited, stack)
                
                # Backtrack: undo changes
                order.pop()
                stack.add(node)
                for u, v in removed_edges:
                    self.DG.add_edge(u, v)
                    if self.DG.in_degree(v) == 1:  # If this node now has incoming edges, remove from stack
                        stack.discard(v)
        
        # Start with nodes that have no incoming edges
        initial_stack = {node for node in self.DG.nodes if self.DG.in_degree(node) == 0}
        backtrack([], set(), initial_stack)

    def var2id(self, variable: Variable):
        """
        Maps a Variable object to a SAT variable index, creating a new mapping if necessary.

        variable: Variable
            The Variable object containing details like type, qubit, col, timeframe, etc.

        Returns:
            int: The SAT variable index for the given Variable.
        """
        if variable not in self.var_bimap.inverse:
            self.var_bimap[self.next_var_index] = variable
            self.next_var_index += 1
        return self.var_bimap.inverse[variable]

    def id2var(self, id: int):
        """
        Extracts the Variable object from the SAT variable index.

        var: int
            The SAT variable index.

        Returns:
            Variable: The corresponding Variable object.
        """
        if id not in self.var_bimap:
            raise ValueError(f"Variable {id} not found in mapping.")
        return self.var_bimap[id]

    def add_variable(self, var_type, **kwargs):
        """
        Creates and maps a Variable object to a SAT variable index.

        Parameters:
            var_type (str): The type of the variable ('state', 'op', 'monitor').
            **kwargs: Additional attributes required to define the variable.
                    - For 'state': qubit, col, timeframe.
                    - For 'op': op_type, qubit, timeframe.
                    - For 'monitor': monitor_type, col, timeframe.

        Returns:
            int: The SAT variable index corresponding to the created Variable.
        """
        # Ensure qubit is a tuple if provided
        qubit = kwargs.get("qubit")
        if isinstance(qubit, list):  # Convert list to tuple
            kwargs["qubit"] = tuple(qubit)

        # Extract attributes based on variable type and validate inputs
        if var_type == "state":
            variable = Variable(
                var_type=var_type,
                axis=kwargs.get("axis"),
                qubit=kwargs.get("qubit"),
                col=kwargs.get("col"),
                timeframe=kwargs.get("timeframe"),
                N=self.qubitNum
            )
        elif var_type == "op":
            variable = Variable(
                var_type=var_type,
                op_type=kwargs.get("op_type"),
                qubit=kwargs.get("qubit"),
                timeframe=kwargs.get("timeframe")
            )
        elif var_type == "monitor":
            variable = Variable(
                var_type=var_type,
                monitor_type=kwargs.get("monitor_type"),
                col=kwargs.get("col"),
                timeframe=kwargs.get("timeframe")
            )
        elif var_type == "aux":
            variable = Variable(
                var_type=var_type,
                aux_type=kwargs.get("aux_type"),
                col=kwargs.get("col"),
                timeframe=kwargs.get("timeframe")
            )
        else:
            raise ValueError(f"Unknown variable type: {var_type}")

        # Map the variable to a SAT variable index
        return self.var2id(variable)
    
    def add_clause(self, clause):
        """
        Adds a single clause directly to the SAT solver and stores it for debugging.
        
        Parameters:
            clause (List[int]): A clause represented as a list of literals (e.g., [1, -2, 3]).
        """
        self.solver.add_clause(clause)
        if not hasattr(self, 'cnf'):  # Ensure `self.cnf` is initialized
            self.cnf = []
        self.cnf.append(clause)    
    
    def print_variables(self):
        """
        Prints all SAT variables and their corresponding Variable objects for debugging.
        Only prints attributes that are not None.
        """
        print("Variables:")
        for var_id, variable in self.var_bimap.items():
            # Filter out None attributes
            variable_attrs = {key: value for key, value in variable.__dict__.items() if value is not None and key != "N"}
            # Format the output
            attrs_str = ", ".join(f"{key}={value}" for key, value in variable_attrs.items())
            print(f"ID: {var_id}, Variable: {attrs_str}")
    
    def print_clauses(self, detail=False):
        """
        Prints all CNF clauses currently stored in the `self.cnf` list.
        
        Parameters:
            detail (bool): If True, prints detailed variable names (e.g., M^t_{ij}).
                        If False, prints simple variable IDs (e.g., x_i).
        """
        if not hasattr(self, 'cnf') or not self.cnf:
            print("No clauses available to print.")
            return

        print("Clauses:")
        for i, clause in enumerate(self.cnf):
            formatted_clause = " ∨ ".join(
                f"{'-' if lit < 0 else ''}{self.var_bimap[abs(lit)].format_flag() if detail and abs(lit) in self.var_bimap else f'x{abs(lit)}'}"
                for lit in clause
            )
            print(f"Clause {i + 1}: ({formatted_clause})")
    
    def equivalence(a, b):
        """
        Generates CNF clauses for logical _equivalence a <=> b.
        
        Parameters:
            a (int): SAT variable for a.
            b (int): SAT variable for b.
        
        Returns:
            List[List[int]]: CNF clauses for a <=> b.
        """
        return [[-a, b], [-b, a]]
    
    def xor_clauses(self, a, b, c):
        """
        Generates CNF clauses for c = a XOR b.
        
        Parameters:
            a (int): SAT variable for a.
            b (int): SAT variable for b.
            c (int): SAT variable for c, representing the result of a XOR b.
        
        Returns:
            List[List[int]]: CNF clauses for c = a XOR b.
        """
        return [
            [-c, -a, -b],  # Clause: ¬c ∨ ¬a ∨ ¬b
            [-c, a, b],  # Clause: ¬c ∨ a ∨ b
            [c, -a, b],  # Clause: c ∨ ¬a ∨ b
            [c, a, -b]     # Clause: c ∨ a ∨ ¬b
        ]
    
    def and_clauses(self, a, b, c):
        """
        Generates CNF clauses for A = B AND C.
        A = B AND C <=> (¬A ∨ B) ∧ (¬A ∨ C) ∧ (A ∨ ¬B ∨ ¬C).

        Parameters:
            a (int): SAT variable for A.
            b (int): SAT variable for B.
            c (int): SAT variable for C.

        Returns:
            List[List[int]]: CNF clauses for A = B AND C.
        """
        return [
            [-a, b],       # ¬A ∨ B
            [-a, c],       # ¬A ∨ C
            [a, -b, -c]    # A ∨ ¬B ∨ ¬C
        ]
    
    def generate_state_matrix(self, timeframe):
        """
        Generates a 2N x M state matrix for the specified timeframe.

        Parameters:
            timeframe (int): The timeframe for which the state matrix is generated.

        Returns:
            List[List[int]]: A 2D list of SAT variable indices representing the state matrix.
        """
        state_matrix = []  # List to store the state matrix
        for row_idx in range(2 * self.qubitNum):
            row = []
            for col_idx in range(self.rotationNum):
                axis = "Z" if row_idx < self.qubitNum else "X"

                var = self.add_variable(
                    var_type="state",
                    axis=axis,
                    qubit=row_idx if row_idx < self.qubitNum else row_idx - self.qubitNum,
                    col=col_idx,
                    timeframe=timeframe
                )
                row.append(var)
            state_matrix.append(row)
        return state_matrix
    
    def buildInitState(self):
        """
        Builds the initial state of the tableau as a SAT formula.
        Encodes the 2N x M tableau into SAT clauses.
        """
        assert(self.current_state == None)
        self.current_state = self.generate_state_matrix(0)  # Generate state matrix for timeframe 0

        for row_idx, row in enumerate(self.tableau):
            for col_idx, value in enumerate(row):
                sat_var = self.current_state[row_idx][col_idx]
                # Add a clause representing the value in the tableau
                if value == 1:
                    self.add_clause([sat_var])  # Clause: var is True
                else:
                    self.add_clause([-sat_var])  # Clause: var is False

    def add_one_hot_constraint(self, vars):
        """
        Adds a global one-hot constraint for a list of operation variables.

        Parameters:
            vars (List[int]): List of SAT variable indices representing all operations at a given timeframe.
        """
        # At least one is true
        self.add_clause(vars)

        # At most one is true
        for i in range(len(vars)):
            for j in range(i + 1, len(vars)):
                self.add_clause([-vars[i], -vars[j]])  # Pairwise negation

    def add_one_depth_constraint(self, operations, timeframe):
        """
        Adds constraints to ensure at most one operation occurs at each qubit per timeframe.

        Constraints:
            1. At least one operation is applied.
            2. S_i and H_i shouldn't both be true for all i.
            3. S_i and CX_jk shouldn't both be true for either i=j or i=k.
            4. H_i and CX_jk shouldn't both be true for either i=j or i=k.
            5. CX_ij and CX_km shouldn't both be true if i or j is included in [k, m].
            6. CX_ij and CX_ji cannot exist at the same time.

        Parameters:
            operations (List[int]): A list of variable IDs representing all operations for the timeframe.
            timeframe (int): The timeframe at which the constraint is applied.
        """

        # Ensure at least one operation is applied
        self.add_clause(operations)

        # Constraint 2: S_i and H_i shouldn't both be true for all i
        for qubit in range(self.qubitNum):
            S_var = self.add_variable(var_type="op", op_type="S", qubit=qubit, timeframe=timeframe)
            H_var = self.add_variable(var_type="op", op_type="H", qubit=qubit, timeframe=timeframe)
            self.add_clause([-S_var, -H_var])  # ¬S_i ∨ ¬H_i

        # Constraints 3 & 4: S_i and H_i shouldn't both be true with CX_jk if i=j or i=k
        for i in range(self.qubitNum):
            for j in range(self.qubitNum):
                if i == j: continue
                S_i = self.add_variable(var_type="op", op_type="S", qubit=i, timeframe=timeframe)
                S_j = self.add_variable(var_type="op", op_type="S", qubit=j, timeframe=timeframe)
                H_i = self.add_variable(var_type="op", op_type="H", qubit=i, timeframe=timeframe)
                H_j = self.add_variable(var_type="op", op_type="H", qubit=j, timeframe=timeframe)
                CX_ij = self.add_variable(var_type="op", op_type="CX", qubit=[i,j], timeframe=timeframe)
                self.add_clause([-S_i, -CX_ij])  # ¬S_i ∨ ¬CX_ij
                self.add_clause([-S_j, -CX_ij])  # ¬S_j ∨ ¬CX_ij
                self.add_clause([-H_i, -CX_ij])  # ¬H_i ∨ ¬CX_ij
                self.add_clause([-H_j, -CX_ij])  # ¬H_j ∨ ¬CX_ij

        cx_conflicts = {i: [] for i in range(self.qubitNum)}  # Track qubits involved in CX
        for i in range(self.qubitNum):
            for j in range(i + 1, self.qubitNum):  # Enforce CX(i, j) only for i < j
                cx_conflicts[i].append([i, j])
                cx_conflicts[i].append([j, i])
                cx_conflicts[j].append([i, j])
                cx_conflicts[j].append([j, i])
        
        # Constraint 5: CX_ij and CX_km shouldn't both be true if they share a or more qubit
        for i in range(self.qubitNum):
            cx_list = list(cx_conflicts[i])  # Get all CX gates involving qubit i
            for idx in range(len(cx_list)):
                CX_ij = self.add_variable(var_type="op", op_type="CX", qubit=cx_list[idx], timeframe=timeframe)
                for jdx in range(idx + 1, len(cx_list)):  # Avoid duplicate checks
                    CX_km = self.add_variable(var_type="op", op_type="CX", qubit=cx_list[jdx], timeframe=timeframe)
                    self.add_clause([-CX_ij, -CX_km])  # ¬CX_ij ∨ ¬CX_km


    def get_property_clauses(self, timeframe, value):
        """
        Generates the clauses needed to encode the property for a given timeframe.

        Property: All columns are "completed" at the specified timeframe.

        Parameters:
            timeframe (int): The current timeframe.
            value (bool): If True, the property represents "all columns are completed."
                        If False, the property represents "at least one column is not completed."

        Returns:
            List[List[int]]: The clauses representing the property.
        """
        # Collect "completed" variables for all columns at the given timeframe
        completed_vars = [
            self.add_variable(var_type="monitor", monitor_type="completed", col=col_idx, timeframe=timeframe)
            for col_idx in range(self.rotationNum)
        ]

        if value:
            # Clauses for "all columns are completed": (C1 ∧ C2 ∧ ... ∧ Cn)
            return [[var] for var in completed_vars]
        else:
            # Clauses for "at least one column is not completed": ¬(C1 ∧ C2 ∧ ... ∧ Cn) = ¬C1 ∨ ¬C2 ∨ ... ∨ ¬Cn
            return [[-var for var in completed_vars]]
    
    def addProperty(self, timeframe):
        """
        Adds the property variable P and its _equivalence to "all columes of the last layer are completed"
        for the specified timeframe.

        Property:
            P <=> (C_i1 ∧ C_i2 ∧ ... ∧ C_in)

        Parameters:
            timeframe (int): The current timeframe.

        Returns:
            int: The SAT variable representing the property P.
        """
        # Create a property variable P for the timeframe
        property_var = self.add_variable(var_type="monitor", monitor_type="property", timeframe=timeframe)

        # Collect the "last layer" columns from the dependency graph
        last_layer = [col_idx for col_idx in range(self.rotationNum) if self.DG.out_degree(col_idx) == 0]
        
        # Collect "completed" variables for only the last layer columns
        last_layer_completed_vars = [
            self.add_variable(var_type="monitor", monitor_type="completed", col=col_idx, timeframe=timeframe)
            for col_idx in last_layer
        ]

        # Add P => (C1 ∧ C2 ∧ ... ∧ Cn)
        for var in last_layer_completed_vars:
            self.add_clause([-property_var, var])  # ¬P ∨ Ci

        # Add (C1 ∧ C2 ∧ ... ∧ Cn) => P
        self.add_clause([property_var] + [-var for var in last_layer_completed_vars])  # P ∨ ¬C1 ∨ ¬C2 ∨ ... ∨ ¬Cn

        return property_var
    
    def assumeProperty(self, timeframe, value):
        """
        Temporarily assumes the property P for the specified timeframe.

        Property:
            P is a variable representing "all columns are completed."

        Parameters:
            timeframe (int): The current timeframe.
            value (bool): True to assume P is true, False to assume P is false.
        """
        # Get the property variable P
        property_var = self.add_variable(var_type="monitor", monitor_type="property", timeframe=timeframe)

        # Assume P is true or false
        self.assumption = [property_var] if value else [-property_var]

    def assertProperty(self, timeframe, value):
        """
        Permanently asserts the property P for the specified timeframe.

        Property:
            P is a variable representing "all columns are completed."

        Parameters:
            timeframe (int): The current timeframe.
            value (bool): True to assert P is true, False to assert P is false.
        """
        # Get the property variable P
        property_var = self.add_variable(var_type="monitor", monitor_type="property", timeframe=timeframe)

        # Assert P is true or false
        if value:
            self.add_clause([property_var])  # Assert P is true
        else:
            self.add_clause([-property_var])  # Assert P is false

    def solve(self, use_assumption=False):
        """
        Solves the SAT formula with an optional flag to use `self.assumption`.

        Parameters:
            use_assumption (bool): If True, uses `self.assumption` for solving.
                                If False, solves without any assumptions.

        Returns:
            bool: True if the formula is satisfiable under the given conditions, False otherwise.
        """
        if use_assumption:
            # Solve the SAT problem with assumptions
            result = self.solver.solve(assumptions=self.assumption)
        else:
            # Solve the SAT problem without assumptions
            result = self.solver.solve()

        # Return whether the formula is satisfiable
        return result

    def save_result(self, timeframe=None):
        """
        Saves the counterexample from the SAT solver if the formula is satisfiable.

        Parameters:
            timeframe (int, optional): The timeframe when the counterexample is found. Default is None.
        """
        # Check if the formula is satisfiable
        if not self.solver.get_model():
            self.saved_result = None  # Clear any previously saved results
            return

        # Extract the model from the solver
        model = self.solver.get_model()
        self.setMaxDepth(timeframe)

        # Decode the model into human-readable format
        counterexample = {}
        for var_id in model:
            if abs(var_id) in self.var_bimap:
                variable = self.id2var(abs(var_id))
                value = var_id > 0  # True if the variable is positive, False if negative
                counterexample[self.var_bimap[abs(var_id)]] = value

        # Save the counterexample for later use
        self.saved_result = {"timeframe": timeframe, "counterexample": counterexample}
    
    def print_result(self, style="detail"):
        """
        Prints the result based on the specified style.

        Parameters:
            style (str): The style of the output. Options are:
                - "cex": Prints the counterexample saved in `self.saved_result`.
                - "detail": Prints detailed tableau resynthesis results.
        """
        # Check if the SAT problem is unsatisfiable
        if not self.solver.get_model():
            print("\nThe SAT problem is unsatisfiable.")
            return

        # Handle different output styles
        if style == "cex":
            self.print_cex()
        elif style == "detail":
            self.print_detail()
        else:
            print(f"Unknown style: {style}. Supported styles are 'cex' and 'detail'.")
    
    def print_cex(self):
        """
        Prints the counterexample saved in `self.saved_result`.
        """
        if not self.saved_result or not self.saved_result.get("counterexample"):
            print("No counterexample found.")
            return

        timeframe = self.saved_result["timeframe"]
        counterexample = self.saved_result["counterexample"]

        print(f"Counterexample Found at Timeframe {timeframe}:")
        for variable, value in counterexample.items():
            # Filter out None configurations in the variable attributes
            filtered_attrs = {k: v for k, v in variable.__dict__.items() if v is not None}
            attrs_str = ", ".join(f"{k}={v}" for k, v in filtered_attrs.items())
            print(f"Variable({attrs_str}): {'True' if value else 'False'}")

    def print_detail(self):
        """
        Prints detailed tableau resynthesis results.
        """
        print("=== Tableau Resynthesis Results ===")
        
        # Iterate through timeframes
        for timeframe in range(self.getMaxDepth() + 1):
            print(f"\nAt timeframe {timeframe}:")

            # Print the state matrix
            state_matrix = self.generate_state_matrix(timeframe)
            for row in state_matrix:
                print("[", end="")
                print(" ".join(str(int(self.solver.get_model()[abs(var) - 1] > 0)) for var in row), end="")
                print("]")

            # Print validity and completion status for each column
            print("\nColumn Properties:")
            for col_idx in range(self.rotationNum):
                valid_var = self.add_variable(var_type="monitor", monitor_type="valid", col=col_idx, timeframe=timeframe)
                completed_var = self.add_variable(var_type="monitor", monitor_type="completed", col=col_idx, timeframe=timeframe)
                valid_status = self.solver.get_model()[abs(valid_var) - 1] > 0
                completed_status = self.solver.get_model()[abs(completed_var) - 1] > 0
                print(f"  Column {col_idx}: valid = {int(valid_status)}, complete = {int(completed_status)}")

            if timeframe == self.getMaxDepth(): break
            
            # Identify the operation applied
            print("\nOperation Applied:")
            operation_found = False
            for i in range(self.qubitNum):
                for op_type in ["S", "H"]:
                    op_var = self.add_variable(var_type="op", op_type=op_type, qubit=i, timeframe=timeframe)
                    if self.solver.get_model()[abs(op_var) - 1] > 0:
                        print(f"  Apply operation: {op_type}_{i}")
                        operation_found = True
                for j in range(self.qubitNum):
                    if i != j:
                        cx_var = self.add_variable(var_type="op", op_type="CX", qubit=[i, j], timeframe=timeframe)
                        if self.solver.get_model()[abs(cx_var) - 1] > 0:
                            print(f"  Apply operation: CX_{i}{j}")
                            operation_found = True

            if not operation_found:
                print("  No operation applied. Something might be wrong")    
      
    def assumeRelease(self):
        """
        Releases all assumptions made in the SAT solver.
        """
        self.assumption = []  # Clear the list of assumptions

    def addMonitor(self, timeframe):
        """
        Adds monitoring constraints for the specified timeframe:
        - Valid Flag: Ensures the column satisfies validity conditions.
        - Completed Flag: Tracks whether the column can be finalized at this timeframe.

        Parameters:
            timeframe (int): The timeframe for which monitoring is added.
        """
        self.addValidFlag(timeframe)  # Add validity constraints
        self.addCompletedFlag(timeframe)  # Add completion constraints

    def addValidFlag(self, timeframe):
        """
        Adds validity constraints for each column in the tableau at the specified timeframe.
        A column is considered valid if:
        - Exactly one Z_i is true (one-hot constraint on Z variables).
        - All X_i are false (ensures X variables for the column are 0).

        This function uses auxiliary variables to modularize and simplify the CNF representation:
        - `alo`: At least one Z variable is true.
        - `amo`: At most one Z variable is true.
        - `ohz`: Combines `alo` and `amo` (i.e., one-hot Z variable).
        - `axf`: All X variables are false.
        - `valid_col`: Indicates that a column satisfies all the above constraints.

        Parameters:
            timeframe (int): The timeframe for which the validity constraints are added.
        """
        # Generate state matrix for the current timeframe
        state_matrix = self.current_state

        for col_idx in range(self.rotationNum):
            # Create a monitor variable for this column
            valid_col = self.add_variable(var_type="monitor", monitor_type="valid", col=col_idx, timeframe=timeframe)

            # Create auxiliary variables
            ohz = self.add_variable(var_type="aux", aux_type="ohz", col=col_idx, timeframe=timeframe)
            alo = self.add_variable(var_type="aux", aux_type="alo", col=col_idx, timeframe=timeframe)
            amo = self.add_variable(var_type="aux", aux_type="amo", col=col_idx, timeframe=timeframe)
            axf = self.add_variable(var_type="aux", aux_type="axf", col=col_idx, timeframe=timeframe)

            # Collect Z_i and X_i for this column
            Z_vars = [state_matrix[i][col_idx] for i in range(self.qubitNum)]  # First N rows are Z_i
            X_vars = [state_matrix[i + self.qubitNum][col_idx] for i in range(self.qubitNum)]  # Last N rows are X_i

            # Define at least one Z_i (alo)
            self.add_clause([-alo] + Z_vars)  # ¬alo ∨ Z_1 ∨ Z_2 ∨ ... ∨ Z_N
            for Z_var in Z_vars:
                self.add_clause([alo, -Z_var])  # alo ∨ ¬Z_i

            # Define "at most one Z_i" (amo) with pairwise constraints and exclude clauses
            for i in range(len(Z_vars)):
                for j in range(i + 1, len(Z_vars)):
                    # ¬amo ∨ ¬Z_i ∨ ¬Z_j
                    self.add_clause([-amo, -Z_vars[i], -Z_vars[j]])

            for exclude_idx in range(len(Z_vars)):
                clause = [amo]  # Start with amo
                for include_idx in range(len(Z_vars)):
                    if include_idx != exclude_idx:
                        clause.append(Z_vars[include_idx])  # Add all other Z_vars
                self.add_clause(clause)

            # Define all X_i are false (axf)
            for X_var in X_vars:
                self.add_clause([-axf, -X_var])  # ¬axf ∨ ¬X_i
            self.add_clause([axf] + X_vars)  # axf ∨ X_1 ∨ X_2 ∨ ... ∨ X_N

            # Define ohz = alo AND amo
            for clause in self.and_clauses(ohz, alo, amo):
                self.add_clause(clause)

            # Define valid_col = ohz AND axf
            for clause in self.and_clauses(valid_col, ohz, axf):
                self.add_clause(clause)

    def addCompletedFlag(self, timeframe):
        """
        Adds constraints for the "completed" flag at the specified timeframe.
        A column is considered "completed" if:
        - At t = 0: It depends only on the validity of the column and its predecessors.
        - At t > 0: It was "completed" in the previous timeframe, OR
                    all its dependencies in the dependency graph are "completed" AND the column is valid.

        Parameters:
            timeframe (int): The timeframe for which the "completed" flag is added.

        Logic:
        - For each column in the tableau, determine its "completed" status by evaluating:
            - Its "completed" status in the previous timeframe (for t > 0).
            - The "completed" status of its predecessors.
            - The validity of the column.
        - Encode these conditions as CNF clauses for the SAT solver.
        """
        for col_idx in range(self.rotationNum):
            # Define the "completed" variable for this column at the current timeframe
            completed_current = self.add_variable(var_type="monitor", monitor_type="completed", col=col_idx, timeframe=timeframe)

            # Define the "valid" variable for this column at the current timeframe
            valid_col = self.add_variable(var_type="monitor", monitor_type="valid", col=col_idx, timeframe=timeframe)

            # Collect dependencies of the current column from the dependency graph
            predecessors_completed = [
                self.add_variable(var_type="monitor", monitor_type="completed", col=dep_idx, timeframe=timeframe)
                for dep_idx in self.DG.predecessors(col_idx)
            ]

            if timeframe > 0:
                # At t > 0, include dependency logic with completed_previous
                completed_previous = self.add_variable(var_type="monitor", monitor_type="completed", col=col_idx, timeframe=timeframe - 1)

                # Add dependency logic
                # completed_current => completed_previous OR (var and valid)
                for var in predecessors_completed:
                    self.add_clause([-completed_current, completed_previous, var])  # completed_current => completed_previous OR var
                self.add_clause([-completed_current, completed_previous, valid_col])   # completed_current => completed_previous OR valid_col

                # completed_previous OR (var and valid) => completed_current
                self.add_clause([completed_current, -completed_previous])
                self.add_clause([completed_current, -valid_col] + [-var for var in predecessors_completed])
            else:
                # At t = 0, no completed_previous; depends only on predecessors_completed and valid_col
                # completed_current => var and valid
                for var in predecessors_completed:
                    self.add_clause([-completed_current, var])  # completed_current => completed_previous OR var
                self.add_clause([-completed_current, valid_col])   # completed_current => completed_previous OR valid_col

                # (var and valid) => completed_current
                self.add_clause([completed_current, -valid_col] + [-var for var in predecessors_completed])

    def addTransition(self, timeframe, Is_obj_depth = True):
        """
        Adds the transition relation for the given timeframe to the SAT formula.

        Parameters:
            timeframe (int): The current timeframe t for which the transition is being added.
            Is_obj_depth (bool): true if set each transition operations to depth = 1
                                else, set it into #gates = 1
        """
        # Step 1: Reuse the current state
        current_state = self.current_state  # Reuse existing current state
        next_state = self.generate_state_matrix(timeframe + 1)  # Generate next state matrix for t+1

        # Step 2: Build operation variables at timeframe t
        operations = []
        for i in range(self.qubitNum):
            # Single-qubit operations
            S_var = self.add_variable(var_type="op", op_type="S", qubit=i, timeframe=timeframe)
            H_var = self.add_variable(var_type="op", op_type="H", qubit=i, timeframe=timeframe)
            operations.append(S_var)
            operations.append(H_var)

            # Two-qubit operations
            for j in range(self.qubitNum):
                if i == j: continue
                CX_var = self.add_variable(var_type="op", op_type="CX", qubit=[i, j], timeframe=timeframe)
                operations.append(CX_var)

        # Step 3: Set depth or gate count constraint
        if Is_obj_depth: self.add_one_depth_constraint(operations, timeframe)
        else: self.add_one_hot_constraint(operations)

        # Step 4: Add constraints for next state variables
        for i in range(self.qubitNum):
            # Rows for Z and X in the current and next state
            Z_i_t = current_state[i]                   # Current Z row for qubit i
            X_i_t = current_state[i + self.qubitNum]   # Current X row for qubit i
            Z_i_t1 = next_state[i]                     # Next Z row for qubit i
            X_i_t1 = next_state[i + self.qubitNum]     # Next X row for qubit i


            # Single-qubit operations
            S_var = self.add_variable(var_type="op", op_type="S", qubit=i, timeframe=timeframe)
            H_var = self.add_variable(var_type="op", op_type="H", qubit=i, timeframe=timeframe)

            # Z_i^{t+1} = Z_i^t ⊕ X_i^t for all columns under S_i
            for col_idx in range(self.rotationNum):
                for clause in self.xor_clauses(Z_i_t[col_idx], X_i_t[col_idx], Z_i_t1[col_idx]):
                    self.add_clause([-S_var] + clause)

            # Z_i^{t+1} = X_i^t, X_i^{t+1} = Z_i^t for all columns under H_i
            for col_idx in range(self.rotationNum):
                self.add_clause([-H_var, Z_i_t1[col_idx], -X_i_t[col_idx]])  # ¬H_var ∨ Z_i_t+1 = X_i_t
                self.add_clause([-H_var, -Z_i_t1[col_idx], X_i_t[col_idx]])  # ¬H_var ∨ Z_i_t+1 = X_i_t
                self.add_clause([-H_var, X_i_t1[col_idx], -Z_i_t[col_idx]])  # ¬H_var ∨ X_i_t+1 = Z_i_t
                self.add_clause([-H_var, -X_i_t1[col_idx], Z_i_t[col_idx]])  # ¬H_var ∨ X_i_t+1 = Z_i_t

            # Two-qubit operations, ctrl: i, target: j
            for j in range(self.qubitNum):
                if i == j: continue
                CX_var = self.add_variable(var_type="op", op_type="CX", qubit=[i, j], timeframe=timeframe)

                # X_j^{t+1} = X_j^t ⊕ X_i^t for all columns
                X_j_t = current_state[j + self.qubitNum]     # Current X row for qubit j
                X_j_t1 = next_state[j + self.qubitNum]       # Next X row for qubit j
                for col_idx in range(self.rotationNum):
                    for clause in self.xor_clauses(X_j_t[col_idx], X_i_t[col_idx], X_j_t1[col_idx]):
                        self.add_clause([-CX_var] + clause)

                # Z_i^{t+1} = Z_i^t ⊕ Z_j^t for all columns
                Z_j_t = current_state[j]  # Current Z row for qubit j
                for col_idx in range(self.rotationNum):
                    for clause in self.xor_clauses(Z_i_t[col_idx], Z_j_t[col_idx], Z_i_t1[col_idx]):
                        self.add_clause([-CX_var] + clause)

            
            # Default: Z_i^{t+1} = Z_i^t, X_i^{t+1} = X_i^t if no operation
            for col_idx in range(self.rotationNum):
                # Collect relevant CX variables where i is the control qubit
                relevant_target_CX_vars = [ # CXs that target at i
                    self.add_variable(var_type="op", op_type="CX", qubit=[_ctrl, i], timeframe=timeframe)
                    for _ctrl in range(self.qubitNum) if _ctrl != i
                ]
                relevant_control_CX_vars = [ # CXs that ctrl at i
                    self.add_variable(var_type="op", op_type="CX", qubit=[i, _tgt], timeframe=timeframe)
                    for _tgt in range(self.qubitNum) if _tgt != i
                ]

                # Z_i^{t+1} = Z_i^t
                self.add_clause([S_var, H_var] + [CX_var for CX_var in relevant_control_CX_vars] + [-Z_i_t1[col_idx], Z_i_t[col_idx]])
                self.add_clause([S_var, H_var] + [CX_var for CX_var in relevant_control_CX_vars] + [Z_i_t1[col_idx], -Z_i_t[col_idx]])

                # X_i^{t+1} = X_i^t
                self.add_clause([H_var] + [CX_var for CX_var in relevant_target_CX_vars] + [-X_i_t1[col_idx], X_i_t[col_idx]])
                self.add_clause([H_var] + [CX_var for CX_var in relevant_target_CX_vars] + [X_i_t1[col_idx], -X_i_t[col_idx]])

            self.current_state = next_state

    def setMaxDepth(self, depth: int):
        """
        set the maximum depth for model checking.
        """
        self.MaxDepth = depth
    
    def getMaxDepth(self):
        """
        Returns the maximum depth for model checking.
        """
        return self.MaxDepth

    def UBMC(self):
        """
        Unbounded Model Checking (UBMC) for tableau resynthesis.
        Handles all timeframes, including t=0, within the loop.
        """
        # Build initial state
        self.buildInitState()
        self.build_rotation_dependency_graph()

        # UBMC loop for all timeframes
        for i in range(self.getMaxDepth()):
            self.addMonitor(i)
            self.addProperty(i)  # Define P <=> (all columns completed at t=i)
            self.assumeProperty(i, True)  # Assume P is true

            # Solve under the assumption that P is true
            if self.solve(use_assumption=True):
                self.save_result(i)  # Save counterexample if SAT
                return

            if i == self.getMaxDepth() - 1:
                # self.assertProperty(i, True)
                self.save_result(i)
                return
            
            # If UNSAT, release assumptions and assert P is false
            self.assumeRelease()  # Clear assumptions
            self.assertProperty(i, False)  # Assert P is false

            # Add transitions for all timeframes except the last one
            
            self.addTransition(i)

    def debug(self):
        """
        For debug
        """
        # Build initial state
        self.buildInitState()
        self.build_rotation_dependency_graph()

        # UBMC loop for all timeframes
        for i in range(self.getMaxDepth()):
            start_time = time.time()  # Start timing for this timeframe
            
            self.addMonitor(i)
            self.addProperty(i)  # Define P <=> (all columns completed at t=i)
            self.assumeProperty(i, True)  # Assume P is true

            # Solve under the assumption that P is true
            if self.solve(use_assumption=True):
                end_time = time.time()  # End timing
                elapsed_time = end_time - start_time
                print(f"✅ Solution found at timeframe {i}. Time taken: {elapsed_time:.4f} seconds.")
                self.save_result(i)  # Save counterexample if SAT
                return

            end_time = time.time()  # End timing
            elapsed_time = end_time - start_time
            print(f"⏳ Timeframe {i} completed, no solution found. Time taken: {elapsed_time:.4f} seconds.")

            if i == self.getMaxDepth() - 1:
                # self.assertProperty(i, True)
                self.save_result(i)
                return
            
            # If UNSAT, release assumptions and assert P is false
            self.assumeRelease()  # Clear assumptions
            self.assertProperty(i, False)  # Assert P is false

            # Add transitions for all timeframes except the last one
            
            self.addTransition(i)


    
            
