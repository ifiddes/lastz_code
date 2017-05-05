import argparse, logging, sys, os, collections, math, random, itertools
import warnings, traceback, multiprocessing
import vcf
import pulp
import scipy.stats, scipy.optimize
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#hard coded region we are looking at (hg19 reference)
region=(120554845,120572645)


class PenaltyTree(object):
    """
    Maintains a tree of penalty terms, so that we can have arbitrarily many
    penalty terms without arbitrarily large constraint expressions.
    
    """
    
    def __init__(self, degree=100):
        """
        Make a new PenaltyTree. degree specifies the maximum number of terms to
        sum together at once. Only one PenaltyTree may be used on a given LP
        problem.
        
        """
        
        # This holds the number of children per node/number of terms to sum at
        # once.
        self.degree = degree
        
        # This holds all our leaf-level terms.
        self.terms = []
        
    def get_variable(self):
        """
        Return a fresh LpVariable with a unique name.
        
        """
        
        # Make the variable
        var = pulp.LpVariable("PenaltyTree_{}".format(get_id()))
        
        # Give the fresh variable to the caller.
        return var
        
    def add_term(self, term):
        """
        Add the given LP expression as a term in the tree.
        
        """
        
        self.terms.append(term)
        
    def set_objective(self, problem):
        """
        Add the sum of all terms as the given LP problem's objective. The
        PenaltyTree must have at least one term.
        
        """
        
        # Algorithm: Go through our leaves, making a variable for the sum of
        # each group of self.degree terms. Then repeat on the sums, making sums
        # of sums, and so on. When we only have one sum, make that the
        # objective.
        
        # This holds the list we're collecting
        collecting = self.terms
        
        # This holds the list of sum variables
        sums = []
        
        while len(collecting) > 1:
            logging.debug("Collecting {} terms in groups of {}".format(len(
                collecting), self.degree))
            
            for i in xrange(0, len(collecting), self.degree):
                # This holds the terms we have collected for this sum
                collected = []
                for j in xrange(0, min(self.degree, len(collecting) - i)):
                    # Grab each term up to our degree that actually exists
                    collected.append(collecting[i + j])
                    
                # This holds the variable we use to represent the sum of the
                # things we have collected
                sum_var = self.get_variable()

                # Constrain this variable to equal the sum of all the collected
                # terms
                problem += sum(collected) == sum_var

                # Add this variable for the next level of collection
                sums.append(sum_var)
               
            # Move up a level in the tree
            collecting = sums
            sums = []
           
        # We have now collected everything down to one term, which is in the
        # collecting list. Use it as our objective function.
        problem += collecting[0]
        
class SequenceGraphLpProblem(object):
    """
    Represents an LP copy number problem. You can attach several models to them,
    constrain them together, and solve the problem.
    
    Internally, contains a pulp LpProblem, and a PenaltyTree.
    
    """
    
    def __init__(self):
        """
        Make a new SequenceGraphLpProblem that we can solve.
        
        """
        
        # We need an actual LpProblem
        self.problem = pulp.LpProblem("copynumber", pulp.LpMinimize)
        
        # We also need a PenaltyTree for organizing penalty terms
        self.penalties = PenaltyTree()
        
    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints to the CopyNumberLpProblem's internal
        LpProblem, and penalties to the model's PenaltyTree.
        
        """
        
        # Make an LP variable for the amount that var_b is above var_a. Note
        # that str() on variables produces their names. Also, we have to make
        # sure that this starts with a letter.
        amount_over = pulp.LpVariable("over_{}".format(get_id()), 0)
            
        # Add the constraint for not being more than that much over
        self.add_constraint(var_b <= var_a + amount_over)
        
        # Make an LP variable for the amount that var_b is below var_a
        amount_under = pulp.LpVariable("under_{}".format(get_id()), 0)
            
        # Add the constraint for not being more than that much under
        self.add_constraint(var_b >= var_a - amount_under)
        
        # Apply an equal penalty in each direction
        self.add_penalty((penalty * amount_over) + (penalty * amount_under))
            
    def add_penalty(self, term):
        """
        Add the given penalty term to the problem's objective.
        
        """
        
        # Just put the term in the PenaltyTree
        self.penalties.add_term(term)
        
    def add_constraint(self, constraint):
        """
        Add the given (exact) constraint to the problem. For approximate
        constraints, use constrain_approximately_equal() instead.
        
        """
        
        # Just add the constraint to the internal problem.
        self.problem += constraint
        
    def solve(self, save=None):
        """
        Solve the LP problem with the best solver we can find. After solving,
        you can get_ploidy() on all the AlleleGroups in Models attached to this
        problem.
        
        If save is specified, it is a filename to which to save the LP problem
        in LP format.
        
        You may only solve a SequenceGraphLpProblem once.
        
        """
        
        # Set up the penalties described by the penalty tree
        logging.info("Setting up penalties")
        self.penalties.set_objective(self.problem)
        
        if save is not None:
            logging.info("Saving problem to {}".format(save))
            self.problem.writeLP(save)
        
        logging.info("Looking for solver")
        # This is a list of solvers to use in reverse order of priority (last is
        # best)
        candidates = [
            pulp.GLPK(),
            # We're not actually using 132 threads here, just 32; the 100 means
            # use deterministic parallelism, working around the segfault issue
            # described here: 
            # <http://list.coin-or.org/pipermail/cbc/2013-March/001044.html>
            pulp.solvers.COIN_CMD(threads=132, msg=1)
        ]
        
        # This is the solver we actually use
        solver = None
        for candidate in candidates:
            if candidate.available():
                logging.info("{} is available".format(
                    candidate.__class__.__name__))
                solver = candidate
            else:
                logging.info("{} is unavailable".format(
                    candidate.__class__.__name__))
        
        if solver is None:
            logging.critical("No solver found!")
            raise Exception("No LP solver available")
        
        logging.info("Solving with {}...".format(solver.__class__.__name__))
        
        # Solve the problem
        status = self.problem.solve(solver)
        logging.info("Solution status: {}".format(pulp.LpStatus[status]))
        
        if len(self.problem.variables()) < 20:
            # It's short enough to look at.
            for var in self.problem.variables():
                logging.debug("\t{} = {}".format(var.name, pulp.value(var)))
                
        # Report the total penalty we got when we solved.
        logging.info("Penalty: {}".format(pulp.value(self.problem.objective)))
        
        if status != pulp.constants.LpStatusOptimal:
            raise Exception("Unable to solve problem optimally.")


class Window(object):
    """Stores the allele fractions of the A and B probes in a specific window
    as well as the LP variables to be solved for this window"""
    def __init__(self, A, B, mid, min_ploidy=0):
        #save the midpoint of the window
        self.mid = mid
        #save the list of A values as a numpy array
        self.A_values = np.array(A)
        #save the list of B values as a numpy array
        self.B_values = np.array(B)
        #save the LP variables that represent the inferred copy number at this window
        self.A = pulp.LpVariable("A", min_ploidy, cat="Integer")
        self.B = pulp.LpVariable("B", min_ploidy, cat="Integer")

    def get_values(self):
        """returns the value lists without removing outliers"""
        return [self.A_values, self.B_values]

    def get_best_estimate(self):
        """Returns a single value for A and B in this window that is a best estimate"""
        A = 10 * np.mean(self.get_values_outliers_removed(A))
        B = 10 * np.mean(self.get_values_outliers_removed(B))
        return A, B

    def get_values_outliers_removed(self):
        """Returns values with outliers removed"""
        return [self.reject_outliers(self.A_values), self.reject_outliers(self.B_values)]

    def reject_outliers(self, data, m = 2.):
    """http://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list"""
    if len(data) > 1:
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d / mdev if mdev else 0.
        return data[s < m]
    else:
        return data

    def get_copy_number(self):
        """returns copy numbers of A and B. LP must be solved"""
        return [pulp.value(self.A), pulp.value(self.B)]


class Model(object):
        """
    Represents a Sequence Graph model of the copy number of two genes. It
    contains a list of Window objects. Approximate-equality constraints
     can be generated tying Window copy numbers to those of their neighbors.
    
    All constraints and penalty terms are automatically added to the
    SequenceGraphLpProblem to which the Model is attached (specified on
    construction).
    
    """
    def __init__(self, problem):
        self.problem = problem
        self.windows = dict()

    def constrain_approximately_equal(self, var_a, var_b, penalty=1):
        """
        Constrain the two LP variables (or constants) var_a and var_b to be
        approximately equal, subject to the given penalty.
        
        Adds the appropriate constraints and penalties to the Model's
        SequenceGraphLpProblem.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        # Just forward the call on to the problem.
        self.problem.constrain_approximately_equal(var_a, var_b, 
            penalty=penalty)
            
    def add_constraint(self, constraint):
        """
        Add the given constraint to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        self.problem.add_constraint(constraint)
    
    def add_penalty(self, penalty): 
        """
        Add the given penalty term to the problem this Model is attached to.
        
        This method is a convenience function, so that you can add constraints
        to a Model when you have only the Model, without poking around inside
        it.
        
        """
        self.problem.add_penalty(penalty)

    def add_window(self, window):
        """Adds a window object to the windows dict, mapped by midpoint"""
        self.windows[window.mid] = window

    def get_windows(self):
        """generator that yields pairs of midpoint, window object"""
        for mid, window in self.windows.iteritems():
            yield mid, window

    def sort_windows(self):
        """Sorts windows by midpoint"""
        self.windows = collections.OrderedDict(sorted(self.windows.items()))

    def build_model(self):
        """Adds measurement aliasing LP variables to this model.
        Returns a dict of LP variables"""
        #holds LP variables by window midpoint
        A_variables = dict()
        B_variables = dict()
        for mid, window in self.get_windows():
            A_value, B_value = window.get_best_estimate()
            A_var = pulp.LpVariable("mid_A_{}".format(mid), 0)
            B_var = pulp.LpVariable("mid_B_{}".format(mid), 0)
            self.add_constraint(A_var == A_value)
            self.add_constraint(B_var == B_value)
            A_variables[mid] = A_var
            B_variables[mid] = B_var
        return A_variables, B_variables

    def add_constraints(self, breakpoint_penalty=1, default_ploidy=2):
        """Constraints adjacent windows together. Breakpoint penalty constraints
        the number of breakpoints allowed. Default ploidy penalizes deviation from
        diploid."""
        for mid, window in self.get_windows():
            #we are in the first window; constrain to default (diploid)
            last_window = None
            if last_window is None:
                self.constrain_approximately_equal(default_ploidy, 
                    window.A, penalty = breakpoint_penalty)
                self.constrain_approximately_equal(default_ploidy, 
                    window.B, penalty = breakpoint_penalty)
            else:
                #constrain this window to be equivalent to last subject to penalty
                self.constrain_approximately_equal(window.A, last_window.A,
                    penalty = breakpoint_penalty)
                self.constrain_approximately_equal(window.B, last_window.B,
                    penalty = breakpoint_penalty)
            last_window = window
        #constrain the very last window to be default as well
        self.constrain_approximately_equal(last_window.A, default_ploidy,
            penalty = breakpoint_penalty)
        self.constrain_approximately_equal(last_window.B, default_ploidy,
            penalty = breakpoint_penalty)

def plot_it(x, A, B, pngpath, samplename):
    """Plot the inferred copy numbers to pngpath
    x = list of mids"""
    x = np.array(x, dtype="int")
    A = np.array(A, dtype="int")
    B = np.array(B, dtype="int")
    fig = plt.figure()
    plt.axis([x[0], x[-1], 0, 5])
    plt.plot(x, A+B, color="blue", label="NOTCH2NL-A")
    plt.plot(x, B, color="red", label="NOTCH2NL-B")
    plt.fill_between(x,A+B, facecolor="blue", alpha=0.7)
    plt.fill_between(x,B, color="red", alpha=0.7)
    r = mpatches.Patch(color="red", label="NOTCH2NL-B", alpha=0.7)
    b = mpatches.Patch(color="blue", label="NOTCH2NL-A", alpha=0.7)
    plt.legend([b, r],["NOTCH2NL-A","NOTCH2NL-B"])
    plt.suptitle("{} NOTCH2NLA/B Copy Number Near Exon 3".format(samplename))
        plt.xlabel("hg19 Position")
    plt.ylabel("Inferred Copy Number")
    plt.savefig(pngpath, format="png")


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--A", type=str, help="A VCF file", required=True)
    parser.add_argument("--B", type=str, help="B VCF file", required=True)
    parser.add_argument("--whitelist", type=str, help="whitelist file with probe weights", default="/hive/users/ifiddes/n2nl_paired_new_suns/whitelist.txt")
    parser.add_argument("--window", type=int, help="window size to use. Default = 5000bp", default=5000)
    parser.add_argument("--step", type=int, help="step size to use. Default = 1000bp", default=1000)
    parser.add_argument("--png", type=str, help="PNG to write out to", required=True)
    parser.add_argument("--name", type=str, help="Patient Name", required=True)
    parser.add_argument("--weight", type=int, help="Breakpoint weight", default=1)
    return parser.parse_args()

def main(args):
    args = parse_args(args)

    #map whitelist positions to paralog and weight (and CHM1 pos)
    wl = [x.split() for x in open(args.whitelist)]
    wl = {x[0]:(x[1],x[2],x[3]) for x in wl}

    #make dict mapping positions to adjusted alelle fractions and positions
    value_dict = dict()
    for v in [args.A, args.B]:
        for record in vcf.Reader(file(v)):
            if str(record.POS) in wl:
                paralog, weight, chm1_pos = wl[str(record.POS)]
                value_dict[record.POS] = paralog, float(weight) * float(record.INFO["ALTFRAC"][0])

    #generate a list of windows that we will calculate the copy number over
    windows = list()
    for x in xrange(region[0], region[1] - args.window, args.step):
        windows.append([x,x + args.window])

    #initialize problem and model
    problem = SequenceGraphLpProblem()
    model = Model(problem)

    #start populating the model with Windows
    for window in windows:
        #find midpoint
        mid = sum(window) / 2
        #make list of values for each paralog
        A = list(); B = list()
        #make list of positions within the current window in the data
        positions = [x for x in value_dict.keys() if x >= window[0] and x <= window[1]]
        #iterate over these and assign the values to A or B
        for position in positions:
            paralog, value = value_dict[position]
            if paralog == "A":
                A.append(value)
            else:
                B.append(value)
        #add these to the model as a new Window
        model.add_window(Window(A, B, mid))

    #time to finish building variables and run model
    model.sort_windows()
    model.add_constraints(args.breakpoint)
    problem.solve()

    x, A, B = list(), list(), list()
    for mid, window in model.get_windows():
        x.append(mid)
        a, b = window.get_copy_number()
        A.append(a); B.append(b)

    plot_it(x, A, B, args.png, args.name)


if __name__ == "__main__" :
    sys.exit(main(sys.argv))