% Clear Command Window
clc;


% Print a blank line
fprintf("\n");


% Print (and enumerate) the information about
% the conditions addressed for the realism no-go theorems
fprintf("First, we need to address the following\n" + ...
        "conditions for the realism no-go theorems:\n");

% Print a blank line
fprintf("\n");

% Print the information and description about
% the Parameter Independence (PI) condition
fprintf(" <strong>1) Parameter Independence (PI):</strong>\n");
fprintf("    * The outcome result of a quantum measurement at\n" + ... 
        "      one location does not depend on the parameter choice of\n" + ...
        "      quantum measurement (e.g., detector orientation) at\n" + ...
        "      a different distant location;\n");

% Print a blank line
fprintf("\n");

% Print the information and description about
% the Outcome Independence (OI) condition
fprintf(" <strong>2) Outcome Independence (OI):</strong>\n");
fprintf("    * The outcome result of a quantum measurement at\n" + ... 
        "      one location does not depend on the outcome result of\n" + ...
        "      another quantum measurement performed at\n" + ...
        "      a different distant location;\n");

% Print a blank line
fprintf("\n");

% Print the information and description about
% the Measurement Independence (MI) condition
fprintf(" <strong>3) Measurement Independence (MI):</strong>\n");
fprintf("    * The configurations of a quantum measurement are\n" + ... 
        "      independently chosen from the quantum properties of\n" + ...
        "      the particles that are being measured;\n");

% Print a blank line
fprintf("\n");

% Print the information and description about
% the Outcome Determinism (OD) condition
fprintf(" <strong>4) Outcome Determinism (OD):</strong>\n");
fprintf("    * The quantum properties of the particles\n" + ...
        "      have well-defined values, independently of\n" + ...
        "      any quantum measurement performed;\n");

% Print a blank line
fprintf("\n");

% Print the information and description about
% the Measurement Non-Contextuality (MNC) condition
fprintf(" <strong>5) Measurement Non-Contextuality (MNC):</strong>\n");
fprintf("    * The value assigned to any quantum property does\n" + ...
        "      not depend on the set of other quantum properties\n" + ...
        "      that are measured simultaneously;\n");


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("***************************************" + ...
        "***************************************");

% Print two blank lines
fprintf("\n\n");


% Print a blank line
fprintf("\n");


% Print (and enumerate) the information about
% the conjunction of conditions for the realism no-go theorems
fprintf("The two no-go theorems addressed in this exercise are:\n");

% Print a blank line
fprintf("\n");

% Print the information and description
% about the conjunction of conditions for
% the Bell - Clauser-Horne-Shimony-Holt (Bell-CHSH) no-go theorem
fprintf(" <strong>1) Bell - Clauser-Horne-Shimony-Holt (Bell-CHSH):</strong>\n");
fprintf("    * <strong>Parameter Independence (PI) ∧ Outcome Independence (OI) ∧\n" + ...
        "      Measurement Independence (MI) ∧ Quantum Mechanics (QM)\n" + ...
        "      => Contradiction!</strong>\n");
fprintf("    * This contradiction occurs because since\n" + ... 
        "      Quantum Mechanics predicts quantum correlations\n" + ...
        "      that violate the Bell-CHSH Inequalities, then one of\n" + ...
        "      the respective conditions should be dropped or violated;\n")

% Print a blank line
fprintf("\n");

% Print the information and description
% about the conjunction of conditions for
% the Kochen - Specker (KS) no-go theorem
fprintf(" <strong>2) Kochen - Specker (KS):</strong>\n");
fprintf("    * <strong>Outcome Determinism (OD) ∧ \n" + ...
        "      Measurement Non-Contextuality (MNC) ∧\n" + ...
        "      Quantum Mechanics (QM) => Contradiction!</strong>\n");
fprintf("    * This contradiction occurs because since\n" + ...
        "      Quantum Mechanics predicts that the outcome results\n" + ...
        "      from the quantum measurements depend of the context of\n" + ...
        "      the other quantum measurements performed, then one of\n" + ...
        "      the respective conditions should be dropped or violated;\n");


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("***************************************" + ...
        "***************************************");

% Print two blank lines
fprintf("\n\n");


% Print a blank line
fprintf("\n");

% Print the headline for the analysis of
% each interpretation of Quantum Mechanics,
% regarding each one of the no-go theorems addressed
fprintf("Now we can analyse how the <strong>De Broglie-Bohm</strong> " + ...
        "<strong>Interpretation</strong>,\n" + ...
        "<strong>Copenhagen Interpretation</strong>, " + ...
        "and <strong>Many World Interpretation</strong>\n" + ...
        "deal with the previous <strong>Bell-CHSH</strong> " + ...
        "and <strong>KS</strong> no-go theorems:\n");

% Print a blank line
fprintf("\n");


% Print the headline for the analysis of
% each interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem
fprintf(" <strong>1) For the Bell-CHSH no-go theorem:</strong>\n");

% Print a blank line
fprintf("\n");

% Print the headline description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem
fprintf("    <strong>i)   De Broglie-Bohm Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics\n" + ...
        "           is explicitly non-local, with the pilot wave\n" + ...
        "           theory or quantum potential used to determine\n" + ...
        "           the trajectory of a particle;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Parameter Independence (PI) condition
fprintf("           <strong>a) Parameter Independence (PI) (Violated):</strong>\n");
fprintf("              * This interpretation is explicitally non-local,\n" + ...
        "                with the pilot wave theory or the quantum potential\n" + ...
        "                influencing the quantum properties of a particle\n" + ...
        "                non-locally, accepting the parameter dependency in\n" + ...
        "                the quantum measurements between distant locations;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Outcome Independence (OI) condition
fprintf("           <strong>b) Outcome Independence (OI) (Violated):</strong>\n");
fprintf("              * This interpretation assumes non-local\n" + ...
        "                quantum correlations between the outcome results\n" + ...
        "                from the quantum measurements being performed\n" + ...
        "                between distant locations;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Measurement Independence (MI) condition
fprintf("           <strong>c) Measurement Independence (MI) (Not Violated):</strong>\n");
fprintf("              * This interpretation assumes non-local\n" + ...
        "                quantum correlations between the outcome results\n" + ...
        "                from the quantum measurements being performed\n" + ...
        "                between distant locations;\n");


% Print two blank lines
fprintf("\n\n");

% Print the headline description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem
fprintf("    <strong>ii)  Copenhagen Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics accepts\n" + ...
        "           quantum non-locality and indeterminism as fundamental\n" + ...
        "           characteristics, thus, considering Quantum Mechanics\n" + ...
        "           as intrinsically probabilistic;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Parameter Independence (PI) condition
fprintf("           <strong>a) Parameter Independence (PI) (Violated):</strong>\n");
fprintf("              * This interpretation accepts implicitly\n" + ...
        "                quantum non-locality, where the outcome results\n" + ...
        "                from the quantum measurements can be affected\n" + ...
        "                instantaneously by the configuration of other\n" + ...
        "                quantum measurements being performed in\n" + ...
        "                another distant location;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Outcome Independence (OI) condition
fprintf("           <strong>b) Outcome Independence (OI) (Violated):</strong>\n");
fprintf("              * This interpretation considers that all\n" + ...
        "                the outcome results of quantum measurements\n" + ...
        "                being performed in distant locations are correlated,\n" + ...
        "                thus, being consistent with the quantum predictions\n" + ...
        "                that violate the Bell Inequalities;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Measurement Independence (MI) condition
fprintf("           <strong>c) Measurement Independence (MI) (Not Violated):</strong>\n");
fprintf("              * This interpretation assumes that all possible\n" + ...
        "                configurations of quantum measurements can be\n" + ...
        "                chosen independently of the quantum properties of\n" + ...
        "                the particles being considered;\n");


% Print two blank lines
fprintf("\n\n");


% Print the headline description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem
fprintf("    <strong>iii) Many World Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics\n" + ...
        "           eliminates the needs for wave-function collapse\n" + ...
        "           and for violating Bell Inequalities, having\n" + ...
        "           quantum non-locality naturally built in;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Parameter Independence (PI) condition
fprintf("           <strong>a) Parameter Independence (PI) (Not Violated):</strong>\n");
fprintf("              * This interpretation considers that all possible\n" + ...
        "                outcome results from a quantum measurement occur in\n" + ...
        "                different and separated branches of the universe or\n" + ...
        "                reality, but the chosen parameters on a location do not\n" + ...
        "                affect the outcome results in another distant location,\n" + ...
        "                within the same branch of the universe or reality;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Outcome Independence (OI) condition
fprintf("           <strong>b) Outcome Independence (OI) (Violated):</strong>\n");
fprintf("              * This interpretation takes into account that\n" + ...
        "                the outcome results of quantum measurements\n" + ...
        "                performed on distant locations are correlated\n" + ...
        "                according to the predictions of Quantum Mechanics,\n" + ...
        "                within each specific branch of the universe or reality;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the Bell-CHSH no-go theorem,
% considering the Measurement Independence (MI) condition
fprintf("           <strong>c) Measurement Independence (MI) (Not Violated):</strong>\n");
fprintf("              * This interpretation assumes that all possible\n" + ...
        "                configurations of quantum measurements are independently\n" + ...
        "                performed on different branches of the universe or reality;\n");



% Print three blank lines
fprintf("\n\n\n");


% Print the headline for the analysis of
% each interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem
fprintf(" <strong>2) For the KS no-go theorem:</strong>\n");

% Print a blank line
fprintf("\n");


% Print the headline description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem
fprintf("    <strong>i)   De Broglie-Bohm Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics it is\n" + ...
        "           a deterministic theory of local hidden variables\n" + ...
        "           and seeks to maintain a notion of classical realism;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Outcome Determinism (OD) condition
fprintf("           <strong>a) Outcome Determinism (OD) (Not Violated):</strong>\n");
fprintf("              * This interpretation maintains the determinism of\n" + ...
        "                the quantum events, since the trajectory of a particle\n" + ...
        "                is guided deterministically by the pilot wave theory or\n" + ...
        "                its respective quantum potential;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the De Broglie-Bohm Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Measurement Non-Contextuality (MNC) condition
fprintf("           <strong>b) Measurement Non-Contextuality (MNC) (Violated):</strong>\n");
fprintf("              * This interpretation should accept the contextuality of\n" + ...
        "                the quantum measurements being performed to be consistent\n" + ...
        "                with Quantum Mechanics, meaning that a value assigned to\n" + ...
        "                some measured quantum property depends on the set of other\n" + ...
        "                quantum measurements being performed simultaneously;\n");


% Print two blank lines
fprintf("\n\n");


% Print the headline description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem
fprintf("    <strong>ii)  Copenhagen Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics broadly\n" + ...
        "           accepts the context of the quantum measurements\n" + ...
        "           performed as a fundamental characteristic;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Outcome Determinism (OD) condition
fprintf("           <strong>a) Outcome Determinism (OD) (Violated):</strong>\n");
fprintf("              * This interpretation rejects the notion of classical\n" + ...
        "                determinism, accepting the intrinsically probabilistic\n" + ...
        "                nature of Quantum Mechanics, where the respective\n" + ...
        "                outcome results do not have pre-defined values\n" + ...
        "                before the quantum measurements occurring;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Copenhagen Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Measurement Non-Contextuality (MNC) condition
fprintf("           <strong>b) Measurement Non-Contextuality (MNC) (Violated):</strong>\n");
fprintf("              * This interpretation accepts that the outcome result of\n" + ...
        "                a quantum measurement depends on its respective context,\n" + ...
        "                i.e., it inherently depends on the set of all the other\n" + ...
        "                quantum measurements performed simultaneously;\n");


% Print two blank lines
fprintf("\n\n");


% Print the headline description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem
fprintf("    <strong>iii) Many World Interpretation:</strong>\n");
fprintf("         * This interpretation of Quantum Mechanics considers\n" + ...
        "           the context of the quantum measurements performed\n" + ...
        "           as a fundamental characteristic of the splitting of\n" + ...
        "           different and parallel branches of the universe or\n" + ...
        "           reality for a specific quantum event;\n");

fprintf("\n");

% Print the information and description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Outcome Determinism (OD) condition
fprintf("           <strong>a) Outcome Determinism (OD) (Not Violated):</strong>\n");
fprintf("              * This interpretation considers that a single\n" + ...
        "                outcome result of the set of possibel ones from\n" + ...
        "                a quantum measurement being performed is deterministic\n" + ...
        "                within each resulting individual and parallel\n" + ...
        "                branch of the universe or reality;\n");

% Print a blank line
fprintf("\n");

% Print the information and description for the analysis of
% the Many World Interpretation of Quantum Mechanics,
% in the perspective of the KS no-go theorem, considering
% the Measurement Non-Contextuality (MNC) condition
fprintf("           <strong>b) Measurement Non-Contextuality (MNC) (Violated):</strong>\n");
fprintf("              * This interpretation takes into account that\n" + ...
        "                each quantum measurement being performed occurs\n" + ...
        "                in a specific context within the respective\n" + ...
        "                resulting branch of the universe or reality;\n");


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("***************************************" + ...
        "***************************************");

% Print two blank lines
fprintf("\n\n");


% Print a blank line
fprintf("\n");


% Print the headline for the conclusions section
fprintf("<strong>Conclusions:</strong>\n");

% Print the information about the conclusions section
fprintf(" * From the two addressed no-go theorems, we can see the KS theorem\n" + ...
        "   as imposing more severe restrictions on classical ontologies,\n" + ...
        "   as it not only rules out local hidden variable theories,\n" + ...
        "   but also any theory that attempts to be non-contextual;\n");
fprintf(" * This notion means that any interpretation of Quantum Mechanics that\n" + ...
        "   seeks to maintain a classical and deterministic realism must accept\n" + ...
        "   the inherent contextuality of quantum measurements, which is\n" + ...
        "   a fundamental constraint that affects possible ontologies.\n");


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("***************************************" + ...
        "***************************************");

% Print two blank lines
fprintf("\n\n");


% Print a blank line
fprintf("\n");