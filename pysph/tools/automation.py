from fnmatch import fnmatch
import glob
import itertools
import json
import os
import shlex
import shutil
import subprocess
import time
import numpy as np


class Task:
    """Basic task to run.  Subclass this to do whatever is needed.

    This class is very similar to luigi's Task class.
    """

    def complete(self):
        """Should return True/False indicating success of task.
        """
        return all([os.path.exists(x) for x in self.output()])

    def output(self):
        """Return list of output paths.
        """
        return []

    def run(self, scheduler):
        """Run the task, using the given scheduler if needed.
        """
        pass

    def requires(self):
        """Return iterable of tasks this task requires.

        It is important that one either return tasks that are idempotent or
        return the same instance as this method is called repeateadly.

        """
        return []


class WrapperTask(Task):
    """A task that wraps other tasks and is done when all its requirements
    are done.
    """
    def complete(self):
        return all(r.complete() for r in self.requires())


class TaskRunner:
    """Run given tasks using the given scheduler.
    """
    def __init__(self, tasks, scheduler):
        """Constructor.

        Parameters
        ----------

        tasks: iterable of `Task` instances.
        scheduler: `pysph.tools.job.Scheduler` instance
        """
        self.tasks = tasks
        self.scheduler = scheduler
        self.todo = []
        for task in tasks:
            self.add_task(task)

    def _show_remaining_tasks(self):
        print("%d tasks pending"%len(self.todo))

    def add_task(self, task):
        if not task.complete():
            self.todo.append(task)
            for req in task.requires():
                self.add_task(req)

    def all_requires_are_done(self, task):
        return all(x.complete() for x in task.requires())

    def run(self, wait=5):
        self._show_remaining_tasks()
        while len(self.todo) > 0:
            to_remove = []
            for i in range(len(self.todo) - 1, -1, -1):
                task = self.todo[i]
                if self.all_requires_are_done(task):
                    to_remove.append(task)
                    print("Running task %s..."%task)
                    task.run(self.scheduler)
            for task in to_remove:
                self.todo.remove(task)

            if len(to_remove) > 0:
                self._show_remaining_tasks()

            if len(self.todo) > 0:
                time.sleep(wait)
        print("Finished!")


class PySPHTask(Task):
    """Convenience class to run a PySPH simulation via an automation
    framework.  The class provides a method to run the simulation and
    also check if the simulation is completed.

    """

    def __init__(self, command, output_dir, job_info=None):
        if isinstance(command, str):
            self.command = shlex.split(command)
        else:
            self.command = command
        self.command += ['-d', output_dir]
        self.output_dir = output_dir
        self.job_info = job_info if job_info is not None else {}
        self.job_proxy = None
        self._copy_proc = None

    ###### Public protocol ###########################################

    def complete(self):
        """Should return True/False indicating success of task.
        """
        job_proxy = self.job_proxy
        if job_proxy is None:
            return self._is_done()
        else:
            return self._copy_output_and_check_status()

    def run(self, scheduler):
        from pysph.tools.jobs import Job
        job = Job(
            command=self.command, output_dir=self.output_dir,
            **self.job_info
        )
        self.job_proxy = scheduler.submit(job)

    def clean(self):
        """Clean out any generated results.

        This completely removes the output directory.

        """
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    ###### Private protocol ###########################################

    def _is_done(self):
        """Returns True if the simulation completed.
        """
        if not os.path.exists(self.output_dir):
            return False
        info_fname = self._get_info_filename()
        if not info_fname or not os.path.exists(info_fname):
            return False
        d = json.load(open(info_fname))
        return d.get('completed')

    def _get_info_filename(self):
        files = glob.glob(os.path.join(self.output_dir, '*.info'))
        if len(files) > 0:
            return files[0]
        else:
            return None

    def _check_if_copy_complete(self):
        proc = self._copy_proc
        if proc is None:
            # Local job so no copy needed.
            return True
        else:
            if proc.poll() is None:
                return False
            else:
                if self.job_proxy is not None:
                    self.job_proxy.clean()
                return True

    def _copy_output_and_check_status(self):
        jp = self.job_proxy
        status = jp.status()
        if status == 'done':
            if self._copy_proc is None:
                self._copy_proc = jp.copy_output('.')
            return self._check_if_copy_complete()
        elif status == 'error':
            cmd = ' '.join(self.command)
            msg = 'On host %s Job %s failed!'%(jp.worker.host, cmd)
            print(msg)
            print(jp.get_stderr())
            raise RuntimeError(msg)
        return False


class Problem(object):
    """This class represents a numerical problem or computational
    problem of interest that needs to be solved.

    The class helps one run a variety of commands (or simulations),
    and then assemble/compare the results from those in the `run`
    method.  This is perhaps easily understood with an example.  Let
    us say one wishes to run the elliptical drop example problem with
    the standard SPH and TVF and compare the results and their
    convergence properties while also keep track of the computational
    time.  To do this one will have to run several simulations, then
    collect and process the results.  This is achieved by subclassing
    this class and implementing the following methods:

     - `get_name(self)`: returns a string of the name of the problem.  All
       results and simulations are collected inside a directory with
       this name.
     - `get_commands(self)`: returns a sequence of (directory_name,
       command_string) pairs.  These are to be exeuted before the
       `run` method is called.
     - `run(self)`: Processes the completed simulations to make plots etc.

    See the `EllipticalDrop` example class below to see a full implementation.

    """
    def __init__(self, simulation_dir, output_dir):
        """Constructor.

        Parameters
        ----------

        simulation_dir : str : directory where simulation output goes.
        output_dir : str : directory where outputs from `run` go.
        """
        self.out_dir = output_dir
        self.sim_dir = simulation_dir

        # Setup the simulation instances in the cases.
        self.cases = None
        self.setup()

    ###### Public protocol ###########################################

    def input_path(self, *args):
        """Given any arguments, relative to the simulation dir, return
        the absolute path.
        """
        return os.path.join(self.sim_dir, self.get_name(), *args)

    def output_path(self, *args):
        """Given any arguments relative to the output_dir return the
        absolute path.
        """
        return os.path.join(self.out_dir, self.get_name(), *args)

    def setup(self):
        """Called by init, so add any initialization here.
        """
        pass

    def get_requires(self):
        """Used by task runners like doit/luigi to run required
        commands.
        """
        base = self.get_name()
        result = []
        for name, cmd, job_info in self.get_commands():
            sim_output_dir = self.input_path(name)
            task = PySPHTask(cmd, sim_output_dir, job_info)
            task_name = '%s.%s'%(base, name)
            result.append((task_name, task))
        return result

    def make_output_dir(self):
        """Convenience to make the output directory if needed.
        """
        base = self.output_path()
        if not os.path.exists(base):
            os.makedirs(base)

    def get_name(self):
        """Return the name of this problem, this name is used as a
        directory for the simulation and the outputs.
        """
        raise NotImplementedError()

    def get_commands(self):
        """Return a sequence of (name, command_string, job_info_dict).

        The name represents the command being run and is used as
        a subdirectory for generated output.

        The command_string is the command that needs to be run.

        The job_info_dict is a dictionary with any additional info to be used
        by the job, these are additional arguments to the
        `pysph.tools.jobs.Job` class. It may be None if nothing special need
        be passed.
        """
        if self.cases is not None:
            return [(x.name, x.command, x.job_info) for x in self.cases]
        else:
            return []

    def get_outputs(self):
        """Get a list of outputs generated by this problem.  By default it
        returns the output directory (as a single element of a list).
        """
        return [self.output_path()]

    def run(self):
        """Run any analysis code for the simulations completed.  This
        is usually run after the simulation commands are completed.
        """
        raise NotImplementedError()

    def clean(self):
        """Cleanup any generated output from the analysis code.  This does not
        clean the output of any nested commands.
        """
        for path in self.get_outputs():
            if os.path.exists(path):
                if os.path.isdir(path):
                    shutil.rmtree(path)
                elif os.path.isfile(path):
                    os.remove(path)


def key_to_option(key):
    """Convert a dictionary key to a valid command line option.  This simply
    replaces underscores with dashes.
    """
    return key.replace('_', '-')

def kwargs_to_command_line(kwargs):
    """Convert a dictionary of keyword arguments to a list of command-line
    options.  If the value of the key is None, no value is passed.

    Examples
    --------

    >>> sorted(kwargs_to_command_line(dict(some_arg=1, something_else=None)))
    ['--some-arg=1', '--something-else']
    """
    cmd_line = []
    for key, value in kwargs.items():
        option = key_to_option(key)
        if value is None:
            arg = "--{option}".format(option=option)
        else:
            arg = "--{option}={value}".format(
                option=option, value=str(value)
            )

        cmd_line.append(arg)
    return cmd_line


def linestyles():
    """Cycles over a set of possible linestyles to use for plotting.
    """
    ls = [dict(color=x[0], linestyle=x[1]) for x in
          itertools.product("kbgr", ["-", "--", "-.", ":"])]
    return itertools.cycle(ls)


class Simulation(object):
    """A convenient class to abstract code for a particular simulation.
    Simulation objects are typically created by ``Problem`` instances in order
    to abstract and simulate repetitive code for a particular simulation.

    For example if one were comparing the elliptical_drop example, one could
    instantiate a Simulation object as follows::

        >>> s = Simlation('outputs/sph', 'pysph run elliptical_drop')

    One can pass any additional command line arguments as follows::

        >>> s = Simlation(
        ...     'outputs/sph', 'pysph run elliptical_drop', timestep=0.005
        ... )
        >>> s.command
        'pysph run elliptical_drop --timestep=0.001'
        >>> s.input_path('results.npz')
        'outputs/sph/results.npz'

    The extra parameters can be used to filter and compare different
    simulations.  One can define additional plot methods for a particular
    subclass and use these to easily plot results for different cases.

    One can also pass any additional parameters to the `pysph.tools.jobs.Job`
    class via the job_info kwarg so as to run the command suitably. For
    example::

        >>> s = Simlation('outputs/sph', 'pysph run elliptical_drop',
        ...               job_info=dict(n_thread=4))

    The object has other methods that are convenient when comparing plots.
    Along with the ``compare_cases``, ``filter_cases`` and ``filter_by_name``
    this is an extremely powerful way to automate and compare results.

    """
    def __init__(self, root, base_command, job_info=None, **kw):
        """Constructor

        Parameters
        ----------

        root: str
            Path to simulation output directory.
        base_command: str
            Base command to run.
        job_info: dict
            Extra arguments to the `pysph.tools.jobs.Job` class.
        **kw: dict
            Additional parameters to pass to command.
        """
        self.root = root
        self.name = os.path.basename(root)
        self.base_command = base_command
        self.job_info = job_info
        self.params = dict(kw)
        self._results = None

    def input_path(self, *args):
        """Given any arguments, relative to the simulation dir, return
        the absolute path.
        """
        return os.path.join(self.root, *args)

    @property
    def command(self):
        return self.base_command + ' ' + self.get_command_line_args()

    @property
    def data(self):
        if self._results is None:
            self._results = np.load(self.input_path('results.npz'))
        return self._results

    def get_labels(self, labels):
        render = self.render_parameter
        if isinstance(labels, str):
            return render(labels)
        else:
            s = [render(x) for x in labels]
            s = [x for x in s if len(x) > 0]
            return r', '.join(s)

    def kwargs_to_command_line(self, kwargs):
        return kwargs_to_command_line(kwargs)

    def get_command_line_args(self):
        return ' '.join(self.kwargs_to_command_line(self.params))

    def render_parameter(self, param):
        """Return string to be used for labels for given parameter.
        """
        if param not in self.params:
            return ''
        value = self.params[param]
        if value is None:
            return r'%s'%param
        else:
            return r'%s=%s'%(param, self.params[param])


def compare_runs(sims, method, labels, exact=None):
    """Given a sequence of Simulation instances, a method name, the labels to
    compare and an optional method name for an exact solution, this calls the
    methods with the appropriate parameters for each simulation.

    Parameters
    ----------

    sims: sequence
        Sequence of `Simulation` objects.
    method: str
        Name of a method on each simulation method to call for plotting.
    labels: sequence
        Sequence of parameters to use as labels for the plot.
    exact: str
        Name of a method that produces an exact solution plot.
    """
    ls = linestyles()
    if exact is not None:
        getattr(sims[0], exact)(**next(ls))
    for s in sims:
        m = getattr(s, method)
        m(label=s.get_labels(labels), **next(ls))

def filter_cases(runs, **params):
    """Given a sequence of simulations and any additional parameters, filter
    out all the cases having exactly those parameters and return a list of
    them.
    """
    def _check_match(run):
        for param, expected in params.items():
            if param not in run.params or run.params[param] != expected:
                return False
        return True

    return list(filter(_check_match, runs))

def filter_by_name(cases, names):
    """Filter a sequence of Simulations by their names.  That is, if the case
    has a name contained in the given `names`, it will be selected.
    """
    return sorted(
        [x for x in cases if x.name in names],
        key=lambda x: names.index(x.name)
    )


############################################################################
# Convenient classes that can be used to easily automate a collection
# of problems.

class SolveProblem(Task):
    """Solves a particular `Problem`. This runs all the commands that the
    problem requires and then runs the problem instance's run method.

    The match argument is a string which when provided helps run only a subset
    of the requirements for the problem.
    """

    def __init__(self, problem, match=''):
        self.problem = problem
        self.match = match
        self._requires = [
            task
            for name, task in self.problem.get_requires()
            if len(match) == 0 or fnmatch(name, match)
        ]

    def output(self):
        return self.problem.get_outputs()

    def run(self, scheduler):
        if len(self.match) == 0:
            self.problem.run()

    def requires(self):
        return self._requires


class RunAll(WrapperTask):
    """Solves a given collection of problems.
    """

    def __init__(self, simulation_dir, output_dir, problem_classes,
                 force=False, match=''):
        self.simulation_dir = simulation_dir
        self.output_dir = output_dir
        self.force = force
        self.match = match
        self.problems = self._make_problems(problem_classes)
        self._requires = [
            SolveProblem(problem=x, match=self.match) for x in self.problems
        ]

    ##### Private protocol  ###############################################

    def _make_problems(self, problem_classes):
        problems = []
        for klass in problem_classes:
            problem = klass(self.simulation_dir, self.output_dir)
            if self.force:
                problem.clean()
            problems.append(problem)
        return problems

    ##### Public protocol  ################################################

    def requires(self):
        return self._requires


class Automator:
    """Main class to automate a collection of problems.

    This processess command line options and runs all tasks with a scheduler
    that is configured using the ``jobs.json`` file if it is present. Here is
    typical usage::

        >>> all_problems = [EllipticalDrop]
        >>> automator = Automator('outputs', 'figures', all_problems)
        >>> automator.run()

    """
    def __init__(self, simulation_dir, output_dir, all_problems):
        """Constructor.

        Parameters
        ----------
        simulation_dir : str
            Root directory to generate simulation results in.
        output_dir: str
            Root directory where outputs will be generated by Problem instances.
        all_problems: sequence of `Problem` classes.
            Sequence of problem classes to automate.
        """
        self.simulation_dir = simulation_dir
        self.output_dir = output_dir
        self.all_problems = all_problems
        self.setup_argparse()

    def setup_argparse(self):
        import argparse
        desc = "Automation script to run simulations."
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument(
            '-f', '--force', action="store_true", default=False, dest='force',
            help='Redo the plots even if they were already made.'
        )
        parser.add_argument(
            '-m', '--match', action="store", type=str, default='', dest='match',
            help="Name of the problem to run (uses fnmatch)"
        )

        all_problem_names = [c.__name__ for c in self.all_problems]
        parser.add_argument(
            '-p', '--problems', action="store", type=str, default="all",
            dest='problems',
            help="specifies problems to run either as a string or a "\
            "set of comma-separated-strings (case-insensitive), "\
            "valid names are %s."%all_problem_names
        )
        self.parser = parser

    def _select_problem_classes(self, problems):
        if problems == 'all':
            return self.all_problems
        else:
            names = [x.strip() for x in problems.lower().split(',')]
            return [cls for cls in self.all_problems
                    if cls.__name__.lower() in names]

    def create_scheduler(self):
        from pysph.tools.jobs import Scheduler

        scheduler = Scheduler(root='.')
        if os.path.exists('jobs.json'):
            scheduler.load('jobs.json')
        else:
            scheduler.add_worker(dict(host='localhost'))
            print("Writing simple jobs.json, please edit it.")
            scheduler.save('jobs.json')
        return scheduler

    def run(self):
        """Start the automation.
        """
        args = self.parser.parse_args()
        problem_classes = self._select_problem_classes(args.problems)
        task = RunAll(
            simulation_dir=self.simulation_dir,
            output_dir=self.output_dir,
            problem_classes=problem_classes,
            force=args.force, match=args.match
        )
        scheduler = self.create_scheduler()
        runner = TaskRunner([task], scheduler)
        runner.run()