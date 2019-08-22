#!/usr/bin/env python
# coding: utf-8


class ScreenEvaluator:
    def __init__(self, model_path, ligand):
        """Initialise the ScreenEvaluator with the screening results file. This populates 
        the self.results object. The Model object is not populated upon initialization, to
        allow rapid inspection of ScreenEvaluator objects.
        
        The self.results object contains the hitlist from the screening experiment as a 
        Pandas dataframe, specifying the Maestro Title, Docking Score, and Decoy status
        of each hit.
        
        Parameters
        ----------
        model_path : str
            Name of the docking output .maegz file, beginning with an entry 
            consisting of a single molecule protein and followed by at least one
            entry consisting of a single molecule ligand
        ligand : str
            Name of the ligand used to load the model. This is passed directly to
            the Model() class method. Options include any string corresponding
            to a ligand title, and 'noligand', which triggers automatic identification
            of the ligand. Alternatively, to bypass model loading, can pass 'nomodel'
        
        """
        # Extract the project table results from the docking output file
        self.input_file = os.path.realpath(model_path)
        self.name = os.path.splitext(os.path.basename(self.input_file))[0]
        self.output_dir = os.path.dirname(self.input_file) + "/"
        
        # Model loading is a slow step, and not required if a hitlist
        # has already been written (usually)
        if not (ligand == 'nomodel'):
            self.model = Model(model_path, ligand)
        
        self.results = ''
        self.metrics = {}
    
    
    def write_hitlist(self):
        """Write the docking hitlist to disk. The hitlist is saved as a .csv
        with the following fields:
        s_m_title : The Maestro title field for each entry.
        r_i_docking_score : Docking score for each hit.
        s_user_RetroScreenRole : type of compound (active/decoy)
        
        This data is loaded into self.results as a Pandas dataframe.
        
        """
        
        with open(self.output_dir + self.name + "-results.csv", 'w') as f:
            extractprotein = subprocess.call(["/home/data/Schrodinger2019_1/utilities/proplister",
                                                "-c",
                                                "-p",
                                                "s_m_title",               
                                                "-p",
                                                "r_i_docking_score",
                                                "-p",
                                                "s_user_RetroScreenRole",
                                                "-p",
                                                "s_user_CompoundSet",
                                                self.input_file],
                                               stdout = f)
                
            
    def compare_score(self, active_score, decoy_score):
        """Function used in AUC calculation to build a score tally based on
        the relationship between an active score and a decoy score. 
        
        Currently, this function is set to assume that lower scores are better.
        Accordingly, the AUC is increased when an active has a lower score 
        than a decoy.
        
        Parameters
        ----------
        active_score : float
            A docking score for an active compound.
        decoy_score : float
            A docking score for a decoy compound.
            
        """
        if active_score < decoy_score:
            return float(1)
        elif active_score > decoy_score:
            return float(0)
        elif active_score == decoy_score:
            return float(0.5)
        
        
    def calc_auc(self, actives_scores, decoys_scores, n_actives, n_decoys):
        """Calculate the ROC AUC for a set of actives and decoys.
        
        This routine scales the AUC number depending on how many of the actives
        and decoys present in the input are present in the score lists.
        The user provides the number of actives and decoys originally used
        as input to the screen, and the number of items in the actives_scores
        and decoys_scores lists are used to determine the scaling of the AUC.
        
        Parameters
        ----------
        actives_scores : list
            List of docking scores for actives
        decoys_scores : list
            List of docking scores for decoys
        n_actives : int
            Number of active compounds used as input to this screen
        n_decoys : int
            Number of decoy compounds used as input to this screen
            
        Returns
        -------
        auc : float
            A scaled AUC that reflects the probability that an active will
            score better than a decoy.
            
        """
        
        n_a_found = float(len(actives_scores))
        n_d_found = float(len(decoys_scores))
        n_a_input = float(n_actives)
        n_d_input = float(n_decoys)

        score_compare_tally = float(0)

        for a in actives_scores:
            for d in decoys_scores:
                score_compare_tally += self.compare_score(a, d)

        auc_raw = float(score_compare_tally / (n_a_found * n_d_found))

        auc_scaled = float((auc_raw * (n_a_found / n_a_input) * (n_d_found / n_d_input)) + ((1 - (n_d_found/n_d_input)) * (n_a_found/n_a_input)))
        
        return auc_scaled

    
    def calc_ef_and_sd(self, threshold, n_actives, n_decoys):
        """Calculate the enrichment factor at the selected threshold.
        The EF is calculated as the ratio between an arbitrarily-selected
        fraction of inactives and the fraction of actives retrieved at that
        fraction of inactives. This definition also corresponds to the 
        gradient of a line constructured from the origin to the ROC-curve value
        coordinate at the selected inactive fraction.
        
        Parameters
        ----------
        threshold : float
            The fraction of inactives at which to determine the EF.
            I.e., EF1% corresponds to a threshold of 0.01.
        n_actives : int
            Number of actives used as input to the screen
        n_decoys : int
            Number of decoys used as input to the screen.
        
        Returns
        -------
        ef : float
            The enrichment factor.
        sd : float
            The analytical error of the enrichment factor
            
        """
        
        for (i, a) in zip(self.roclist_x, self.roclist_y):
            if i > threshold:
                ef = float(a/i)
                break
        
        s = (a*np.log(a)) / (i*np.log(i))
        sd = (((a*(1-a))/float(n_actives)) + ((s**2)*i*(1-i))/n_decoys) / (i**2)
        
        return (ef, sd)
    
    
    def calc_ci95(self, variance):
        """Calculates a 95% confidence interval for any metric, with the 
        assumption that the distribution is normal. 
        
        Caveat from the original paper:
        "This will break down for systems with very large errors (generally
        due to having very few actives or decoys) ... When large errors do
        arise the exact magnitude of the error is generally less relevant than
        the fact that the error is large"
        
        Parameters
        ----------
        variance : float
            The variance of any metric
        
        Returns
        -------
        ci95 : float
            The 95% confidence interval. Adding and subtracting this value 
            from the mean value of the metric will give the values that define
            the borders of 95% confidence for that value.
            
        """
        ci95 = 1.96 * variance**0.5
        return ci95
    
    
    def precision_recall(self, hitlist):
        """Calculate the precision-recall profile and associated summary statistics
        (F1 score, Matthews correlation quality, Youden's J, PR-AUC).
        
        The summary statistics are defined for every point along the precision-recall
        curve, so will report the maximal value.
        
        Parameters
        ----------
        hitlist : list
            A list containing 'active' and 'decoy' entries. Can be generated from the self.results.role
            Pandas sheet column.
        
        """
        
        self.precision = []
        self.recall = []
        
        self.f1 = []
        self.mcc = []
        self.youdenj = []
        
        
        # Move the threshold through the hitlist
        for i in range(1, len(hitlist)):
            toplist = hitlist[0:i]
            bottomlist = hitlist[i:]
            
            tp = 0
            fp = 0
            tn = 0
            fn = 0
        
            # At each step, calculate the precision and recall
            for x in toplist:
                if x == 'active':
                    tp += 1
                elif x == 'decoy':
                    fp += 1
            for x in bottomlist:
                if x == 'decoy':
                    tn += 1
                elif x == 'active':
                    fn += 1
            
            #print len(toplist), len(bottomlist), "TP: ", tp, "FP: ", fp, "TN: ", tn, "FN: ", fn
            
            precision = float(tp) / float(tp + fp)
            recall = float(tp) / float(tp + fn)
            
            self.precision.append(precision)
            self.recall.append(recall)
            
            # Calculate the summary statistics
            f1 = np.reciprocal(((np.reciprocal(recall)) + (np.reciprocal(precision)))/2)
            
            self.f1.append(f1)
            
            mcc = float((tp * tn) - (fp * fn)) / float(math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
            self.mcc.append(mcc)
            
            youdenj = (float(tp)/float(tp+fn)) + (float(tn)/float(tn+fp)) - 1
            self.youdenj.append(youdenj)
            
        
        # Calculate PRAUC
        self.metrics['prauc'] = metrics.auc(self.recall, self.precision)
        
        # Store the summary statistics
        self.metrics['f1'] = max(self.f1)
        self.metrics['mcc'] = max(self.mcc)
        self.metrics['youden_j'] = max(self.youdenj)
        
        # Store the score cutoff associated with the maximum Youden's J index
        self.metrics['youden_j_score'] = self.results['score'][self.youdenj.index(max(self.youdenj))]
        
    
    def calc_enrichment(self, n_actives, n_decoys):
        """Calculate enrichment metrics for this screen. The metrics are stored in the 
        self.metrics dictionary object.
        
        Parameters
        ----------
        n_actives : int
            The number of actives in the input library used to perform the screen.
        n_decoys : int
            The number of decoys in the input library used to perform the screen.
            
        """
        # Load hitlist
        try:
            self.results = pd.read_csv(self.output_dir + self.name + "-results.csv", names = ['title', 'score', 'role', 'compoundset'], skiprows = [0,1])
        except IOError:
            print("File not found: " + self.output_dir + self.name + "-results.csv")
            return
        
        # Remove duplicate items from the hitlist.
        self.results.drop_duplicates(subset='title', keep="first", inplace=True)
        
        # Extract score lists for actives and decoys
        self.actives_scores = self.results[self.results["role"] == 'active'].score.tolist()
        self.decoys_scores = self.results[self.results["role"] == 'decoy'].score.tolist()
        self.roclist = self.results.role.tolist()
        
        # Prepare ROC lists
        x = [0]
        y = [0]
        d_count = 0
        a_count = 0
        for i in self.roclist:
            if i == 'decoy':
                d_count += 1
            elif i == 'active':
                a_count += 1
            x.append(d_count)
            y.append(a_count)

        #Scale the coordinates by the total numbers
        self.roclist_x = np.asarray(x, dtype=float)/float(d_count)
        self.roclist_y = np.asarray(y, dtype=float)/float(a_count)
        
        # Populate the self.enrichment object with the metrics
        self.metrics['AUC'] = self.calc_auc(self.actives_scores, self.decoys_scores, n_actives, n_decoys)
        
        q1 = self.metrics['AUC'] / float(2-self.metrics['AUC'])
        q2 = (2*(self.metrics['AUC'])**2) / (1 + self.metrics['AUC'])
        self.metrics['AUC_SD'] = (self.metrics['AUC']*(1-self.metrics['AUC']) + ((n_actives-1)*(q1-(self.metrics['AUC'])**2)) + ((n_decoys-1)*(q2-(self.metrics['AUC'])**2))) / float(n_actives * n_decoys)
        
        self.metrics['EF1'], self.metrics['EF1_SD'] = self.calc_ef_and_sd(0.01, n_actives, n_decoys)
        self.metrics['EF2'], self.metrics['EF2_SD'] = self.calc_ef_and_sd(0.02, n_actives, n_decoys)
        self.metrics['EF5'], self.metrics['EF5_SD'] = self.calc_ef_and_sd(0.05, n_actives, n_decoys)
        self.metrics['EF10'], self.metrics['EF10_SD'] = self.calc_ef_and_sd(0.10, n_actives, n_decoys)
        
        self.metrics['AUC_CI95'] = self.calc_ci95(self.metrics['AUC_SD'])
        self.metrics['EF1_CI95'] = self.calc_ci95(self.metrics['EF1_SD'])
        self.metrics['EF2_CI95'] = self.calc_ci95(self.metrics['EF2_SD'])
        self.metrics['EF5_CI95'] = self.calc_ci95(self.metrics['EF5_SD'])
        self.metrics['EF10_CI95'] = self.calc_ci95(self.metrics['EF10_SD'])
    
    
    def plot_roc(self, color):
        """Plot the ROC curve and save to a file named <self.name>_roc.png.
        
        """

        x = [0]
        y = [0]
        d_count = 0
        a_count = 0
        for i in self.roclist:
            if i == 'decoy':
                d_count += 1
            elif i == 'active':
                a_count += 1
            x.append(d_count)
            y.append(a_count)

        #Scale the coordinates by the total numbers
        x_scale = np.asarray(x, dtype=float)/float(d_count)
        y_scale = np.asarray(y, dtype=float)/float(a_count)

        plt.figure(figsize=(3,3))

        sns.set(context="notebook", style="ticks")

        plt.plot([0,1], [0,1], color=sns.xkcd_rgb["grey"], linestyle='--')
        plt.plot(x_scale, y_scale, sns.xkcd_rgb[color])

        plt.ylim(0,1)
        plt.ylabel("Actives")

        plt.xlim(0,1)
        plt.xlabel("Inactives")

        plt.savefig(self.output_dir + self.name + "_roc.png", dpi=300, format="png", bbox_inches="tight", transparent=True)

    def plot_roc_marked(self, color, watch_compound, watch_role):
        """Plot the ROC curve and add text labels to compounds
        specified.
        
        Parameters
        ----------
        color : str
            The colour to plot the ROC curve in.
        watch_compound : str
            Text label to match against the "compoundset" property for
            selecting the marked compounds.
        watch_role : str
            Text label to match against the "role" property for selecting
            the marked compounds.
        """

        x = [0]
        y = [0]
        marked_y = []
        marked_x = []
        marked_labels = []

        d_count = 0
        a_count = 0

        for i in range(0, len(self.results)):
            if (self.results.iloc[i]['role'] == 'decoy'):
                d_count += 1
            elif (self.results.iloc[i]['role'] == 'active'):
                a_count += 1
            if (self.results.iloc[i]['compoundset'] == watch_compound) and (self.results.iloc[i]['role'] == watch_role):
                marked_y.append(a_count)
                marked_x.append(d_count)
                marked_labels.append(self.results.iloc[i]['title'])

            x.append(d_count)
            y.append(a_count)


        #Scale the coordinates by the total numbers
        x_scale = np.asarray(x, dtype=float)/float(d_count)
        y_scale = np.asarray(y, dtype=float)/float(a_count)
        marked_x_scale = np.asarray(marked_x, dtype=float)/float(d_count)
        marked_y_scale = np.asarray(marked_y, dtype=float)/float(a_count)

        plt.figure(figsize=(6,6))

        sns.set(context="notebook", style="ticks")

        plt.plot([0,1], [0,1], color=sns.xkcd_rgb["grey"], linestyle='--')
        plt.plot(x_scale, y_scale, sns.xkcd_rgb[color])

        plt.scatter(marked_x_scale, marked_y_scale, color='k', marker='x', s=50)

        texts = []
        align = 'right'
        for index,(i,j) in enumerate(zip(marked_x_scale, marked_y_scale)):
            if (align == 'right'):
                texts.append(plt.text(i+.02, j-0.01, marked_labels[index]))
                align = 'left' 
            elif (align == 'left'):
                texts.append(plt.text(i-0.02*len(marked_labels[index]), j-0.01, marked_labels[index]))
                align = 'right'



        plt.ylim(0,1)
        plt.ylabel("Actives")

        plt.xlim(0,1)
        plt.xlabel("Inactives")

        adjust_text(texts)

        plt.savefig(self.output_dir + self.name + "_roc_" + watch_role + "_" + watch_compound + ".png", dpi=300, format="png", bbox_inches="tight", transparent=True)
        plt.close()
        
        
        
    def plot_pr(self, color):
        """Plot the Precision-Recall curve.
        
        """
        
        plt.figure(figsize=(7,7))

        sns.set(context="notebook", style="ticks")
        
        f_scores = np.linspace(0.2, 0.8, num=4)
        for f_score in f_scores:
            x = np.linspace(0.01, 1)
            y = f_score * x / (2 * x - f_score)
            l, = plt.plot(x[y >= 0], y[y>=0], color='gray', alpha=0.2)
            plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

        plt.plot(self.recall, self.precision, sns.xkcd_rgb[color])

        plt.ylim(0,1)
        plt.ylabel("Precision")

        plt.xlim(0,1)
        plt.xlabel("Recall")

        plt.savefig(self.output_dir + self.name + "_pr.png", dpi=300, format="png", bbox_inches="tight", transparent=True)
        plt.close()
        


