/** \mainpage
 *
 *  \image html theta-medium.png
 * 
 * This is the doxygen documentation of %theta, a program for template-based statistical modeling and inference, focussing on problems
 * in high-energy physics.
 * 
 * This doxygen documentation covers the C++ code and the (closely related) configuration file format of theta.
 * 
 * For most users, the recommended way of using theta is via the higher-level, python-based
 * user interface <em>theta-auto</em>, which is documented <a href="http://theta-framework.org/">here</a>.
 *
 * \section main_intro Documentation overview
 *
 * <ol>
 *   <li>\subpage intro Introduction guides you through an example and explains the typical workflow with %theta</li>
 *   <li>\subpage cmd_interface "Command line interface" described the command line tools of %theta, namely the \c theta
 *     program and the \c merge program.</li>
 *   <li>\subpage extend "Extending theta" describes how to extend %theta using the plugin system</li>
 *   <li>\subpage test "Testing" describes how %theta is tested.</li>
 * </ol>
 *
 * \section license License
 *
 * %theta is licensed under the <a href="http://www.gnu.org/copyleft/gpl.html">GPL</a>.
 *
 * \section ref References
 * %theta includes algorithms based on:
 * <ol>
 *   <li><em>Matsumoto, Makoto and Nishimura, Takuji:</em> <a href="http://doi.acm.org/10.1145/272991.272995">"Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator"</a>,
 *        ACM Trans. Model. Comput. Simul. 1, 1998</li>
 *   <li><em>Pierre L'Ecuyer:</em> "Maximally Equidistributed Combined Tausworthe Generators", Math. Comp. 65, 1996</li>
 *   <li><em>Pierre L'Ecuyer:</em> "Tables of Maximally Equidistributed Combined LFSR Generators", Math. Comp. 68, 1999</li>
 *   <li><em>George Marsaglia and Wai Wan Tsang:</em> "The Ziggurat Method for Generating Random Variables", Journal of Statistical Software 8, 2000</li>
 *   <li><em>A. Gelman, G. O. Roberts, and W. R. Gilks:</em> "Efficient Metropolis Jumping Rules", Bayesian Statistics 5, 1996</li>
 * </ol>
 *
 * %theta contains software from:
 * <ul>
 * <li><a href="http://www.hyperrealm.com/libconfig/libconfig.html">libconfig</a> A
 *      well-written, well-documented C/C++ library for processing configuration files with a very simple and elegant API.</li>
 * <li><a href="http://www.chokkan.org/software/liblbfgs/">liblbfgs</a> An implementation of
 *     the Limited-memory Broyden-Fletcher-Goldfarb-Shanno minimization algorithm</li>
 * </ul>
 * These libraries are included in the %theta subversion repository.
 *
 * Furthermore, some parts of numerical algorithms have been copied from  <a href="http://www.gnu.org/software/gsl/">GNU Scientific Library (GSL)</a>.
 */

/**
 * \page intro Introduction
 *
 * This page discusses a concrete example where you get an overview over how %theta works
 * from the point of view of a user. In the second section, some internals of %theta are
 * explained which are good to know even if you do not plan to extend %theta.
 *
 * \section first_example First example
 *
 * Suppose you search for a new particle. After a sophsticated event selection,
 * you have events containing candidates of your new particle.
 * For each of these events, you can reconstruct the mass. From your Monte-Carlo simulation,
 * you conclude that your signal has a distribution in this reconstructed mass in the
 * form of a gaussian with mean \f$ 1000\,\mathrm{GeV}/c^2 \f$ and width \f$ 250\,\mathrm{GeV}/c^2 \f$, whereas your background
 * is expected to be flat in the region from 500 to \f$ 1500\,\mathrm{GeV}/c^2 \f$ which should be used to further constrain
 * your background.
 *
 * You do the studies mainly at one fixed integrated luminosity L. From a background fit
 * to a sideband you expect that you can constrain your background poisson mean in the signal
 * region to \f$ 1600 \pm 200 \f$ events. The model of your signal allows for a large variety of signal
 * cross sections; the standard model predicts no "signal".
 *
 * For this introduction, we consider the creation of likelihood-ratio test statistic distribution for the
 * "background only" null hypothesis in order to find out the critical region for the hypothesis test attempting
 * to reject the null hypothesis at some confidence level.
 *
 * Analysis with %theta always consists of several steps, namely:
 *<ol>
 * <li>Model definition</li>
 * <li>Configuration of the statistical methods to apply</li>
 * <li>Configuration of the run</li>
 * <li>Executing the %theta main program (optional: more than once)</li>
 * <li>(optional:) merging the output produced by different runs of %theta</li>
 * <li>analyzing the output by making plots</li>
 *</ol>
 * All but the last point are well supported by %theta and explained below in more detail.
 *
 * \subsection model_def Model definition
 *
 * The first step to do in any analysis with %theta is to translate you model (like the one above) into a %theta
 * configuration.
 *
 * Now, have a look at the <tt>examples/gaussoverflat.cfg</tt> configuration file. The syntax is actually quite simple
 * but might require some time to get used to. For now, only two
 * recurring terms are introduced: "setting" is any statement of the form "parameter = value;" and
 * "setting group" which is a set of settings enclosed in curly braces, e.g., in
 * <pre>
 * mass = {
 *   range = (500.0, 1500.0);
 *   nbins = 200;
 * };
 * </pre>
 * the right hand side of the "mass" setting is a setting group containing a
 * "range" setting (which has a list as type) and "nbins"
 * setting (an integer type). See the libconfig reference linked in the \ref ack section
 * for a detailed description of the configuration file syntax. Note that any right hand side
 * where a setting group is expected, it is also allowed to write a string of the form "@&lt;path&gt;"
 * where &lt;path&gt; is the (dot-separated) path in the configuration file. In this case, the path is resolved and the value
 * of this setting is used instead. This means, that the following two settings are equivalent:
 * \code
 * //1.
 * double_value = 1.0;
 * setting_group_1 = {
 *   some_value = 1.0;
 * };
 *
 * //2.
 * double_value = "@setting_group_1.some_value";
 * setting_group_1 = {
 *   some_value = 1.0;
 * };
 * \endcode
 *
 * This is especially useful for large models which require a lot of definitions and would be very hard to understand without
 * meaningful links.
 *
 * At the top of the configuration file, the parameters and observables you want to use are defined:
 * there is one observable "mass" with the range [500, 1500] and 200 bins. Note that %theta does
 * not care at all about units and that observables are <em>always</em> binned.
 * The parameters of this model are defined next: "s" is  the (poisson) mean number of signal events after your selection, "b" is the mean number of
 * background events. Only parameter names are defined at this point. The values they take for pseudo data
 * generation are defined later.
 *
 * Next, the model "gaussoverflat" is defined, where the expectation for the "mass" observable is specified. As discussed above,
 * it is a linear combination of s signal events which are gaussian and b background events which are flat. This linear combination
 * of different components is expressed as different setting groups where you specify \c coefficient-function and \c histogram for each component.
 *
 * After the observable specification, there is a special setting named "parameter-distribution". This defined
 * the distribution of the model parameters and will be used
 * <ul>
 * <li>as additional term in the likelihood function (Bayesianically, these are priors)</li>
 * <li>when throwing pseudo experiments as distribution for the model parameters</li>
 * </ul>
 *
 * \subsection conf_stat Configuration of the statistical methods to apply
 *
 * After having defined the model, the next thing to take care of is the statistical method you want to apply.
 *
 * In this case, this is done by the "hypotest" settings group.
 *
 * This settings group defines a statistical method (also called "producer", as it produces
 * data in the output %database) of type "deltanll_hypotest". This producer (which is documented \link deltanll_hypotest here \endlink)
 * expects two more settings: "signal-plus-background-distribution" and "background-only-distribution".
 * They are both setting groups which specify special distributions (which replace the distribution in the model) to apply to
 * get these two model varaints from the original model. In this case, the "background-only" model is given
 * by the fixing the signal parameter s to zero.
 *
 * Whenever a producer modeule runs, it will be provided with a model and data. The \ref deltanll_hypotest producer will
 * <ol>
 * <li> Construct the likelihood function of the model, given the data</li>
 * <li> Minimize the negative logarithm of the likelihood function with the "signal-plus-background" parameter values fixed</li>
 * <li> Minimize the negative log-likelihood function with the "background-only" parameter values fixed</li>
 * <li> save the two values of the negative-log-likelihood in the result table</li>
 * </ol>
 * For more details, refer to the documentation of \ref deltanll_hypotest.
 *
 * \subsection conf_run Configuration of the run
 *
 * Having configured the model and a statistical method, we have to glue them together. In this case, we want to make
 * many pseudo experiments in order to determine the distribution of the likelihood-ration test-statistics.
 *
 * This is done by the "main" settings group. It defines a theta::run which
 * throws random pseudo data from a model and calls a list of producers. The results will be written
 * to what is configured via the "output_database" setting. Currently, you can shoose between \link sqlite_database \endlink
 * and \link rootfile_database \endlink.
 *
 * If you issue the %theta command, it will read the configuration file and build an instance of \link theta::Run \endlink using the "main" setting group.
 * Then, the run method of \link theta::Run theta::Run \endlink method is called which:
 * <ol>
 *  <li>Uses the "data_source" to get the (pseudo-)data in form of a \link theta::Data Data \endlink instance (which is basically one Histogram for each observable).
 *    For the type "model_source", this means generating random values for model parameters and poisson data around the model prediction. But you can also
 *    use a "histo_source" if you have actual data.</li>
 *  <li>The data from the previous step is passed to each of the producers, along with the configured "model".</li>
 *  <li>What the producers do with the \link theta::Model Model \endlink and \link theta::Data Data\endlink is up to them. Usually, they construct the likelihood function
 *    and do something with it. In this case, a producer of the type \link deltanll_hypotest deltanll_hypotest\endlink is run which minimizes the likelihood function twice, using
 *    different parameter distributions.</li>
 *  <li>The results from the producers just run is appended to the "products" table in the configured "output_database". In case of an error, no row is added to
 *    the "products" table. Instead, an entry in the "log" table is done.</li>
 * </ol>
 *
 * These steps are repeated "n-events" times.
 *
 * \subsection running_theta Executing theta
 *
 * To actually do all the work, you have to call %theta with the configuration file name as argument
 * <pre>
 * bin/theta examples/gaussoverflat.cfg
 * </pre>
 *
 * This will produce a SQLite output file called gaussoverflat.db in the directory you called %theta.
 * You can open the file with
 * <pre>
 *  sqlite3 gaussoverflat.db
 * </pre>
 *
 * You will be presented a prompt in which you can enter SQL commands. To see which tables and columns are defined,
 * enter
 * <pre>
 * .schema
 * </pre>
 *
 * You will see that there are three tables:
 * <ul>
 *  <li>'log' which contains log messages of all errors and warnings during the run</li>
 *  <li>'rndinfo' which contains the random number seed used for the run</li>
 *  <li>'products'. This is the most important table which contains one line per pseudo experiment where the results of all producers is saved</li>
 * </ul>
 *
 * While the first two tables always have the same schema, the columns of the last table is defined by the producers which run.
 * In this case, the value of s and b as used in the data source are saved as well as the result of the deltanll_hypotest producer,
 * i.e., \c hypotest__nll_b and \c hypotest__nll_sb. Note that the column names are always built from the 'name' of the producer in the
 * configuration file, followed by two underscores, followed by a producer-specific column name.
 * Using the producer's name as columns name makes it possible to run the same type of producer more than once in one run of %theta.
 * To find out about the columns a producer generates, read the according class documentation (the C++ class name is the name given in the 'type' setting).
 *
 * \subsection analyzing Analyzing the output
 *
 * The output of a run of %theta is saved as SQLite %database to the output file configured in the run
 * settings group. What to do with these parameters very much depends on the use case.
 *
 * %theta only supports very rudimentary root histogram production through SQL queries. For this example, execute
 * <pre>
 * bin/histos examples/plot-gaussoverflat.cfg
 * </pre>
 *
 * This will produce the file \c gaussoverflat.root which contains the histograms as specified in the plotting
 * configuration file. Open this file with root; you will find the following histograms:
 * <ul>
 *  <li>\c log_likelihood_ratio contains the negative logarithm of the likelihood ratio. It is the difference of the \c hypotest__nll_b  and
 *      \c hypotest__nll_sb columns of the \c products table.</li>
 *  <li>\c diced_s, the value for \c s used to create pseudo data. In this case, it is always 200.</li>
 *  <li>\c diced_b, the value for \c b used to create pseudo data. In this case, it is a normal distribution around 1000 with width 200.</li>
 * </ul>
 *
 * \section plugins Available Plugins
 *
 * Everywhere in the configuration file where there is a "type=..." setting, the plugin system is used
 * to construct an object the corresponding type. The string given in the "type=..." setting is the C++ class
 * name used. The configuration is documented in the respective class documentation.
 *
 * To find out about the available plugins, look at the <a href="hierarchy.html">class hierarchy</a>
 * (use "Classes", then "Class Hierarchy" in the above menu). The
 * abstract base classes the plugins are derived from are:
 * <ul>
 * <li>\link theta::HistogramFunction HistogramFunction\endlink: used in the "histogram=..."-setting in the observables
 *      specification of a model.</li>
 * <li>\link theta::Function Function \endlink: used as coefficients of the components of an observable or as prior
 *         specification in a model. The only core plugin available is \link mult mult \endlink, which multiplies a list of parameters.</li>
 * <li>\link theta::Minimizer Minimizer\endlink: used by some producers such as maximum likelihood, profile likelihood methods.</li>
 * <li>\link theta::Distribution Distribution \endlink: used in model constraints and as priors in a statistical method. Used in "main.model.parameter-distribution".
 *     Also, many producers and DataSources allow to use a setting "override-parameter-distribution". </li>
 * <li>\link theta::Producer Producer \endlink: statistical method called by theta::Run which are provided Data and Model and write
 *      results to the theta::ProductsTable. Used in "main.producers".</li>
 * <li>\link theta::DataSource DataSource \endlink: provides the data/pseudo data on which to run the producers. Used in "main.data_source".</li>
 * <li>\link theta::Database Database \endlink: saves the products produced by the producers on disk. Used in "main.output_database".</li>
 * </ul>
 */
 
 
 /**
  * \page extend Extending theta
  *
  * One important concept introduced earlier is the plugin architecture of %theta.
  * As a rule of thumb, every setting group in a %theta configuration file will give rise to one instance
  * of a C++ class during runtime of %theta. More specifically, <b>any setting group containing a <em>type="<typename>"</em>
  * setting is used to construct a C++ object of class &lt;typename&gt; via the plugin system.</b>
  *
  * To define and use your own plugin, you have to:
  * <ol>
  * <li>Define a new class which is derived from a plugin-class and implement all its pure virtual methods and a constructor
  *     taking a \link theta::Configuration Configuration \endlink object as the only argument.</li>
  * <li>In a .cpp-file, call the REGISTER_PLUGIN(yourclass) macro</li>
  * <li>Make sure to compile and link this definition to a shared-object file.</li>
  * <li>In the configuration file, make sure to load the shared-object file as plugin.
  *     You can now use the plugin defined as any other %theta component via
  *     a setting group containing type="yourclass";</li>
  * </ol>
  * For examples of plugins, you can have a look at the \c plugins/ directory, which contains the core plugins of %theta compiled as \c lib/core-plugins.so.
  * The easiest way to get started is to create the new ".cpp" file in the \c plugins/ directory. This way, it will be automatically compiled
  * into \c lib/core-plugins.so and you do not have to care about 3. and 4.
  *
  * \section extend_example Example: Function plugin
  *
  * To have a simple but complete example, consider you want to use a 1/sqrt(s) prior for a Bayesian method.
  * In the corresponding producer, it is possible to pass this prior as the setting "additional-nll-term". This is a
  * setting group which specifies a \link theta::Function Function \endlink, so the plugin required here is a class derived from Function.
  *
  * The first thing to do is to write down the documentation, of how this plugin will be configured. Note that
  * the setting group you speficy has to provide enough information to completely specify the instance to create:
  * the constructor will <em>only</em> have accedd to the \link theta::Configuration Configuration \endlink object which essentially
  * contains the settings from the configuration file and the defined observables and parameters, but not more.
  *
  * The setting group should look like this:
  * \code
  *   some_name = {
  *     type = "nl_one_over_sqrt";
  *     parameter = "s";
  *   };
  * \endcode
  * The "type" setting is always there. This is required by the plugin system to delegate the construction of such an object
  * to the right plugin. The rest of the settings are up to the author of the plugin. Note that, depending on the type of the plugin,
  * base classes also might require some settings.
  *
  * For the implementation, it is considered good style to use the plugin "type" name as filename and to
  * use this name as C++ class name. You should also create a separate header file which contains the class declaration and
  * the documentation (this  is required for the documentation generation). For sake of simplicity, I will omit this
  * complication.
  *
  * The code required for this plugin is:
  * \code
  *
  * using namespace theta;
  * using namespace theta::plugin;
  *
  * class nl_one_over_sqrt: public Function{
  * private:
  *    ParId pid;
  * public:
  *    nl_one_over_sqrt(const Configuration & cfg): pid(cfg.vm->getParId(cfg.setting["parameter"])){
  *    }
  *
  *   virtual double operator()(const ParValues & values) const{
  *      return 0.5 * values.get(pid);
  *   }
  *
  *  };
  *
  * REGISTER_PLUGIN(nl_one_over_sqrt)
  * \endcode
  *
  * Let's go through this step-by-step.
  *
  * The first lines just include the namespaces \c theta and \c theta::plugin in the lookup in order to save
  * some typing.
  *
  * The class \c nl_one_over_sqrt is defined as derived class of \link theta::Function Function \endlink. Each class used in the
  * plugin system must be derived from one of the classes suited for plugins (see section "Available Plugin" on the \subpage intro Introduction
  * page.
  *
  * The \c nl_one_over_sqrt class defines a private data member of type \c ParId. \c ParId instances are used to
  * manage parameter identities: each parameter defined in the configuration file corresponds to a certain
  * value of \c ParId. \c ParId instances can be copied, assigned and compared, but this is all what you can
  * do with them. Any additional information about the parameter (currently only its name)
  * is managed with an instance of \link theta::VarIdManager \endlink.
  *
  * The public constructor must be defined to have a \code const Configuration & \endcode as the only parameter.
  * \c cfg.setting represents the configuration file setting group the instance is constructed from.
  * This can be used to access the settings in this setting group through
  * the \c operator[]. See the documentation of \link theta::SettingWrapper SettingWrapper \endlink and the
  * libconfig documentation for details.
  * The string configured read into \c parameter_name is then used to get the corresponding \c ParId value.
  * This has to be done with an instance of \ref theta::VarIdManager. The currently relevant instance
  * is available through the Configuration instance, \c cfg.vm.
  *
  * Note that the code is a bit complicated by the fact that there is no default constructor for the class \c ParId. This
  * is a design choise meant to prevent access to invalid (i.e., non-existing) parameters. The <em>only</em> way
  * to get instances of \c ParId and \c ObsId is though an instance of \c VarIdManager, usually \c cfg.vm.
  *
  * Before going on, it is also important to see what's <em>not</em> there:
  * <ul>
  *   <li>There is no code which handles the case that there is no "parameter" setting in \c cfg.setting. There does not have
  *      to be as the call \c cfg.setting["parameter"] will throw an exception which will be caught in the %theta
  *      main program and reported to the user.</li>
  *   <li>There is no explicit type casting of the configuration file setting. This is done implicitely via the overloaded
  *     casting operators of \link theta::SettingWrapper SettingWrapper \endlink (this concept is actually
  *     stolen form libconfig). In case the 
  *     required type does not matched the one in the configuration file, an exception will be thrown, which will be caught
  *     in %theta.</li>
  *   <li>There is no code to handle the case that the user has mis-spelled the parameter name in the configuration.
  *      In this case ... yes, you guessed right: an exception will be thrown by VerIdManager::getParId.</li>
  * </ul>
  *
  * So these two lines of construction code are enough to have robust code with (implicitely provided) proper error handling.
  *
  * Next, we have a look at the \c operator(). This implements the purely virtual
  * specification theta::Function::operator() and is what is called if the function value is requested.
  * From the argument \c values, we can request the value of \c pid and return the result. However, this is
  * not quite correct yet: the function does not provide proper error handling in case the parameter value is
  * smaller than zero. Therefore, this code should be modified such that the method body is:
  * \code
  *   double val = values.get(pid);
  *   if(val < 0.0) throw MathException("nl_one_over_sqrt: negative argument");
  *   return 0.5 * val;
  * \endcode
  *
  * The macro REGISTER_PLUGIN registers the new class at the plugin system. You have to call this macro exactly
  * once for every plugin you want to register. Do this in a .cpp file, never in the header file.
  *
  * One piece is still missing: as the documented in theta::Function::par_ids, this variable
  * should contain all the parameters this function depends on.
  *
  * The completed example (including the split of .hpp and .cpp) is implemented in
  * \c plugins/nl_one_over_sqrt.hpp and \c plugins/nl_one_over_sqrt.cpp.
  */
 
 
 /** \page cmd_interface Command line interface
  *
  * %theta comes with two programs meant to be run from the command line: %theta and merge.
  *
  * \section cmd_theta theta
  *
  * The %theta command expects the configuration file(s) to use as argument(s). theta will use the
  * setting group "main" to construct and configure a \link theta::Run \endlink instance.
  * Additionally to this setting group, the special setting group "options" is parsed.
  * This setting group must contain the setting \c plugin_files,
  * which is a list of shared-object filenames to load as plugins. If relative paths are used, they
  * are resolved from where theta is executed.
  *
  * %theta parses the configuration file, constructs the \link theta::Run Run \endlink object from the specified setting group
  * and invoke its \link theta::Run::run run method \endlink.
  *
  * %theta will output the current progress if running on a tty, unless disabled via the \c -q (or \c --quiet) option.
  *
  * Using the \c --nowarn option will disable warnings about unused configuration file settings. Only
  * use this if you are sure that the warnings can be savely ignored (in case of a problem, <em>always</em>
  * reproduce it without using this option first).
  *
  * If you send the \c SIGINT signal to %theta (e.g., by hitting ctrl+C on a terminal running %theta),
  * it will exit gracefully as soon as the current toy experiment is processed. This feature is useful for
  * interactive use if the whole run takes too long but you still want to be able to
  * analyze the output produced so far. It can also be used in a batch job script
  * which can send this signal to theta just before the job reaches the maximum time by which it would be 
  * killed by the batch system.
  *
  * \section cmd_merge merge
  *
  * The \c merge program is used to merge result databases from different runs of the %theta command into a single one.
  * The \c merge program will only work on <em>compatible</em> databases, i.e., databases which contain the same tables
  * and where tables contain the same columns. Usually, you should only attempt to merge result databases which were produced
  * using the same configuration file.
  *
  * \c merge re-maps the runid entry in the input files such that the output file contains unique runid entries.
  *
  * \c merge will refuse to merge result databases which use the same random number generator seed.
  *
  * A typical invocation is :
  * \code
  *  merge --outfile=merged.db result1.db result2.db result3.db
  * \endcode
  * which merges the files result{1,2,3}.db into the output file merged.db. If the output file exists,
  * it will be overwritten.
  *
  * Another possibility is to merge all files matching \c *.db within one directory via
  * \code
  *  merge --outfile=merged.db --in-dir=results
  * \endcode
  *
  * The only other supported option is \c -v or \c --verbose which increases the verbosity of merge.
  */

/** \page testing Testing
 *
 * There are two types of tests in %theta: unit tests which focus on testing the functionality of
 * a single class or method and statistical test cases which compare results obtained from %theta
 * with results analytically obtainable. For the unit tests, look at the source files in the \c test source
 * directory. To run the unit tests, execute
 * <pre>
 * make
 * source setenv.sh
 * make run-test
 * </pre>
 * from the %theta directory.
 *
 * In the following sections, only the statistical tests will be discussed. The source code of these
 * tests can be found in the \c test/test-stat directory.
 *
 * \section testing_counting-nobkg Counting experiment without background
 *
 * As first test case, consider a model with poisson signal around a true mean \f$ \Theta \ge 0 \f$. The outcome of an experiment is completely
 * described by the number of observed events, \f$ n \f$. The likelihood function is given by
 * \f[
 *   L(\Theta | n) = \frac{\Theta^n e^{-\Theta}}{n!}
 * \f]
 * which is maximized by the choice \f$ \Theta = n \f$.
 *
 * A model template is defined in \c test-stat/counting-nobkg.cfg.tpl; the tests will be done for \f$ \Theta=5.0 \f$ (called "low n case" below)
 * and the for \f$ \Theta=10000.0 \f$ (called "asymptotic case" below). Some methods are not expected to have the
 * desired properties for these models. However, this is not what is to be tested here; rather, the implementation
 * of the algorithms is the interest of these tests.
 *
 * \subsection testing_counting-nobkg_mle Maximum likelihood method
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-nobkg-mle.py</tt><br>
 * <em>Relevant plugin classes:</em> \link mle mle \endlink, \link root_minuit root_minuit \endlink
 *
 * 500 pseudo experiments are performed by throwing a Poisson random number around \f$ \Theta \f$ which is passed
 * as number of observed events \f$ n \f$ to the mle implementation.
 *
 * The maximum likelihood method should always estimate \f$ \hat\Theta = n \f$. A small deviation due to numerical evaluation
 * is tolerated. The relative deviation must not exceed \f$ 10^{-4} \f$. In the asymptotic case, also the error estimate
 * \f$ \hat\sigma_{\Theta} \f$ from the minimizer is checked: \f$ \hat\sigma_{\Theta}^2 = \hat\Theta \f$ should hold with
 * a maximum relative error of \f$ 10^{-4} \f$.
 *
 * \subsection testing_counting-nobkg_mcmc_quantiles Bayesian Intervals with Markov-Chain Monte-Carlo 
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-nobkg-mcmc_quant.py</tt><br>
 * <em>Relevant plugin classes:</em> \link mcmc_quantiles mcmc_quantiles \endlink
 *
 * With a flat prior on the signal mean \f$ \Theta \f$, the posterior in
 * case of \f$ n \f$ observed events, \f$ \pi(\Theta | n)\f$, is the same as the likelihood function given above
 * (it is already properly normalized):
 * \f[
 *   \pi(\Theta | n) = \frac{\Theta^n e^{-\Theta}}{n!}.
 * \f]
 *
 * Given an estimate for the q-quantile, \f$ \hat{\Theta}_q \f$, the true value of \f$ q \f$ is given by
 * \f[
 *  q = \int_0^{\hat{\Theta}_q} d\Theta\,\pi(\Theta | n) = \frac{\gamma(n+1, \hat{\Theta}_q)}{\Gamma(n+1)}
 * \f]
 * where \f$ \Gamma \f$ is the gamma function and \f$ \gamma \f$ the incomplete lower gamma function.
 *
 * The tested quantiles correspond to the median, the one-sigma central interval and a 95% upper limit: 0.5, 0.16, 0.84, 0.95.
 * The Markov chain length is 10,000.
 *
 * The test calculates the error on \f$ q \f$ (<em>not</em> on \f$ \hat{\Theta}_q \f$) for 50 pseudo experiments for each of
 * the 4 quantiles. If the deviation is larger than 0.015 (absolute), then the pseudo-experiment
 * is marked as "estimate too low" or "estimate too high"
 * depending whether the requested value of \f$ q \f$ was lower or high than the estimated value. (A pseudo-experiment
 * can have both marks as there are 4 quantiles which are tested.)
 *
 * If there are more than 15 pseudo experiments in either of the categories,
 * the test is considered failed, otherwise, it is considered passed.
 * As additional diagnostics, the number of pseudo experiments marked "too low"
 * and "too high" are printed. In a bias-free method, these
 * two values should be similar, i.e., the value estimated should be sometimes too
 * low and sometimes too high, and not prefer
 * one side. However, no automatic diagnostics is run, as these numbers are small.
 *
 * \subsection testing_counting-nobkg_deltanll "Delta-log-likelihood" confidence intervals
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-nobkg-deltanll_intervals.py</tt><br>
 * <em>Relevant plugin classes:</em> \link deltanll_intervals deltanll_intervals \endlink, \link root_minuit root_minuit \endlink
 *
 * Given the number of observed events, \f$ n \f$, the logarithm of the likelihood ratio between
 * a free value \f$ \Theta \f$ and the value which maximizes the likelihood function (i.e., \f$ \Theta = n \f$) is
 * \f[
 *    \log LR(\Theta | n) = n \log\frac n \Theta - n + \Theta.
 * \f]
 * The minimum is at \f$ \Theta = n \f$ where this logarithm becomes 0.
 *
 * For a given confidence level \f$ c \f$, the interval construction finds the two values of \f$ \Theta \f$,
 * \f$ l_c \f$ and \f$ u_c \f$ with \f$ l_c \le u_c \f$ for which the difference of the log-likelihood to its minimum,
 * \f$ \Delta NLL \f$ corresponds (via the asymptotiv property of the likelihood ratio acording to Wilks' Theorem)
 * to the requested confidence level \f$ c \f$ via the equation
 * \f[
 *   \Delta NLL(c) = \left( \mathrm{erf}^{-1}(c) \right)^2
 * \f]
 * where \f$ \mathrm{erf} \f$ is the error function.
 *
 * For fixed \f$ n \f$ and \f$ c \f$, the values for \f$ \Theta \f$ which yields this difference in negative-log-likelihood is found
 * <ol>
 *  <li>by solving the equation \f$ \log LR(\Theta | n) = \Delta NLL(c)\f$ directly (which has to be done numerically), and</li>
 *  <li>by the generic \link deltanll_intervals deltanll_intervals \endlink method implemented in %theta.</li>
 * </ol>
 *
 * Technically, point 1. is accomplished by an independent python function which expects \f$ n \f$ and \f$ c \f$ as
 * input and applies the bisection method to find the lower and upper value of \f$ \Theta \f$ for which the equation of 1. holds.
 *
 * The requested confidence levels tested correspond to one-sigma and two-sigma interals and the maximum likelihood value, i.e.,
 * confidence levels 0.6827, 0.9545 and 0. The outcome if the method are estimates for the lower and upper values of
 * the intervals for \f$ |Theta \f$. Below, they will be called \f$ l_{1\sigma}, u_{1\sigma}, l_{2\sigma}, u_{2\sigma}, l_{0} \f$.
 *
 * In both cases, the maximum tolerated deviation in the lower and upper interval borders is \f$ 10^{-4}\Theta \f$ .
 *
 * In the asymptotic case, it is also checked that the "one-sigma" and "two-sigma" intervals are where expected in this regime, i.e.,
 * \f{eqnarray*}
 *   (l_{1\sigma} - l_0)^2 &= l_0 \\
 *   (u_{1\sigma} - l_0)^2 &= l_0 \\
 *   (l_{2\sigma} - l_0)^2 &= 4\cdot l_0 \\
 *   (u_{2\sigma} - l_0)^2 &= 4\cdot l_0 \\
 * \f}
 * with a maximum relative difference of \f$ 2\cdot 10^{-2} \f$ , where "relative difference" means dividing the difference of the
 * left and right hand side of these equations by the right hand side (note: the tolerance is chosen relatively large to avoid
 * reporting false positives. This test is mainly to ensure that there is no common bug in the python function for point 1.
 * and the \link deltanll_intervals deltanll_intervals \endlink plugin when calculating \f$ \Delta NLL(c) \f$ ).
 *
 * \section testing_counting-fixedbkg Counting experiment with fixed background mean
 *
 * As second test case, consider a counting experiment with poisson signal with true mean \f$ \Theta \ge 0 \f$ and
 * background mean \f$ \mu \ge 0 \f$ which is assumed to be known. This example is mainly useful to test the correct
 * treatment of a fixed parameter and the error estimates on these parameters.
 *
 * The tests are very similar to the previous section. As there, the values 5.0 and 10000.0 are tested for \f$ \Theta \f$.
 * The values for \f$ \mu \f$ used are 3.0 and 2000.0. For each test case, all four possible combinations of these
 * parameter values are tested.
 *
 * \subsection testing_counting-fixedbkg_mle Maximum likelihood method
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-fixedbkg-mle.py</tt><br>
 * <em>Relevant plugin classes:</em> \link mle mle \endlink, \link root_minuit root_minuit \endlink
 *
 * 500 pseudo experiments are performed by throwing a Poisson random number around \f$ \Theta + \mu \f$ which is passed
 * as number of observed events \f$ n \f$ to the mle implementation.
 *
 * The maximum likelihood method should estimate \f$ \hat\Theta = n - \mu \f$ for \f$ n > \mu \f$ and \f$ \hat\Theta = 0 \f$ otherwise.
 * The maximum tolerated numerical deviation from this result is set to \f$ 10^{-4} \cdot (\Theta + \mu) \f$.
 *
 * In the asymptotic case (\f$ \Theta = 10000 \f$), the error estimate for the minimizer, \f$ \hat{\sigma}_{\Theta} \f$ is checked
 * to fulfill \f$ \hat{\sigma}_\Theta^2 = \Theta + \mu \f$ with a maximum relative error of \f$ 10^{-4} \f$.
 *
 * \subsection testing_counting-fixedbkg_mcmc_quantiles Bayesian Intervals with Markov-Chain Monte-Carlo
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-fixedbkg-mcmc_quant.py</tt><br>
 * <em>Relevant plugin classes:</em> \link mcmc_quantiles mcmc_quantiles \endlink
 *
 * With a flat prior on the signal mean \f$ \Theta \f$ and fixed \f$ \mu \f$, the posterior in case of \f$ n \f$ observed events,
 * \f$ \pi(\Theta | n)\f$, is
 * \f[
 *   \pi(\Theta | n) = \frac{1}{\Gamma(n+1, \mu)} (\Theta + \mu)^n e^{- \Theta - \mu}
 * \f]
 * where \f$ \Gamma \f$ is the "upper" incomplete Gamma function.
 *
 * Given an estimate for the q-quantile, \f$ \hat{\Theta}_q \f$, the true value of \f$ q \f$ is given by
 * \f[
 *  q = \int_0^{\hat{\Theta}_q} d\Theta\,\pi(\Theta | n) = 1 - \frac{\Gamma(n+1, \hat{\Theta}_q + \mu )}{\Gamma(n+1, \mu)}.
 * \f]
 *
 * The tested quantiles correspond to the median, the one-sigma central interval and a 95% upper limit: 0.5, 0.16, 0.84, 0.95.
 * The Markov chain length is 10,000.
 *
 * The test calculates the error on \f$ q \f$ for 50 pseudo experiments for each of
 * the 4 quantiles. If the deviation is larger than 0.015 (absolute), then the pseudo-experiment
 * is marked as "estimate too low" or "estimate too high"
 * depending whether the requested value of \f$ q \f$ was lower or high than the estimated value. (A pseudo-experiment
 * can have both marks as there are 4 quantiles which are tested.)
 *
 * If there are more than 15 pseudo experiments in either of the categories,
 * the test is considered failed, otherwise, it is considered passed.
 * As additional diagnostics, the numbers of pseudo experiments marked "too low"
 * and "too high" are printed.
 *
 * \subsection testing_counting-fixedbkg_deltanll "Delta-log-likelihood" confidence intervals
 *
 * <em>Test script:</em> <tt>test/test-stat/counting-fixedbkg-deltanll_intervals.py</tt><br>
 * <em>Relevant plugin classes:</em> \link deltanll_intervals deltanll_intervals \endlink, \link root_minuit root_minuit \endlink
 *
 * As in case of no background, the likelihood ratio of the likelihood function at \f$ \Theta \f$ and at its minimum is calculated.
 * Now, there are two cases to treat seperately:
 * <ol>
 *  <li> \f$ n >= \mu \f$. In this case, the maximum likelihood estimate for \f$ \Theta \f$ is given by \f$ \hat{\Theta} = n - \mu \f$ </li>
 *  <li> \f$ n < \mu \f$. In this case, the maximum likelihood estimate for \f$ \Theta \f$ is given by \f$ \hat{\Theta} = 0 \f$ </li>
 * </ol>
 *
 * In the first case, the logarithm of the likelihood ratio is given by
 * \f[
 *    \log LR(\Theta | n) = n \log\frac{n}{\Theta + \mu} - n + \Theta + \mu.
 * \f]
 * The minimum is at \f$ \Theta = n - \mu \f$ where this logarithm becomes 0.
 *
 * In the second case, the logarithm of the likelihood ratio is given by
 * \f[
 *  \log LR(\Theta | n) = n \log\frac{\mu}{\Theta + \mu} + \Theta
 * \f]
 *
 *
 * For a given confidence level \f$ c \f$, the interval construction finds the two values of \f$ \Theta \f$,
 * \f$ l_c \f$ and \f$ u_c \f$ with \f$ l_c \le u_c \f$ for which the difference of the log-likelihood to its minimum,
 * \f$ \Delta NLL \f$ corresponds (via the asymptotiv property of the likelihood ratio acording to Wilks' Theorem)
 * to the requested confidence level \f$ c \f$ via the equation
 * \f[
 *   \Delta NLL(c) = \left( \mathrm{erf}^{-1}(c) \right)^2
 * \f]
 * where \f$ \mathrm{erf} \f$ is the error function.
 *
 * For fixed \f$ n \f$ and \f$ c \f$, the values for \f$ \Theta \f$ which yields this difference in negative-log-likelihood is found
 * <ol>
 *  <li>by solving the equation \f$ \log LR(\Theta | n) = \Delta NLL(c)\f$ directly (which has to be done numerically), and</li>
 *  <li>by the generic \link deltanll_intervals deltanll_intervals \endlink method implemented in %theta.</li>
 * </ol>
 *
 * Technically, point 1. is accomplished by an independent python function which expects \f$ n \f$ and \f$ c \f$ as
 * input and applies the bisection method to find the lower and upper value of \f$ \Theta \f$ for which the equation of 1. holds.
 *
 * The requested confidence levels tested correspond to one-sigma and two-sigma interals and the maximum likelihood value, i.e.,
 * confidence levels 0.6827, 0.9545 and 0. The outcome if the method are estimates for the lower and upper values of
 * the intervals for \f$ |Theta \f$. Below, they will be called \f$ l_{1\sigma}, u_{1\sigma}, l_{2\sigma}, u_{2\sigma}, l_{0} \f$.
 *
 * In both cases, the maximum tolerated deviation in the lower and upper interval borders is \f$ 10^{-4}\Theta \f$.
 *
 *
 * \section testing_onoff The "On/Off" and related problems
 *
 * The tests of the following tests use the "on/off" problem and the "gaussian mean" problem are
 * discussed in detail in <a href="http://arxiv.org/abs/physics/0702156">arXiv:0702156</a>, called
 * "Reference" in the following subsections.
 *
 * \subsection testing_onoff_pl On/Off problem with profile likelihood
 *
 * <em>Test script:</em> <tt>test/test-stat/onoff_pl.py</tt><br>
 * <em>Relevant plugin classes:</em> \link deltanll_hypotest deltanll_hypotest \endlink
 *
 * Relevant for this test is in particular Section 5, equation (20) of the Reference.
 *
 * In this test, the Z values given there in Table 1 of the Reference for the profile likelihood method are compared to the output
 * of the deltanll_hypotest producer. The rounded value of the Z value must be identical to the values given in the paper
 * to consider this test passed (if rounding to 2 decimal digits).
 *
 * \subsection testing_gaussmean_pl Gaussian mean problem with profile likelihood
 *
 * <em>Test script:</em> <tt>test/test-stat/gaussmean_pl.py</tt><br>
 * <em>Relevant plugin classes:</em> \link deltanll_hypotest deltanll_hypotest \endlink
 *
 * See equation (21) from Section 5 of the Reference.
 *
 * Z values obtained from deltanll_hypotest are compared to the Z values in Table 1 of the Reference. The test
 * is considered passed iff all rounded Z values agree (if rounding to 2 decimal digits).
 */

/** \brief Common namespace for %theta
 */
namespace theta{}
