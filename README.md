# spinodal-decomposition-probability_ANN
This repository provides MATLAB code for calculating descriptors from alloy compositions and using a trained artificial neural network (ANN) to predict the probability of spinodal decomposition. The workflow enables computational screening of candidate compositions before experimental validation.

**Before you start**
Keep the code and Excel files in the same folder and set it as the MATLAB current folder. The supplied workbook illustrates the Al–Co–Cr–Fe–Ni system. To study another system, prepare your own composition workbook and update the filenames, sheet names, and row ranges in the relevant scripts. MATLAB R2019a or later is required; normalization and training use Deep Learning Toolbox.

**Generate Candidate Compositions**
For a five-element system, run Generate_5fraction.m after setting the composition limits and step size for the selected elements. Use a suitable generator for a different number of elements. Prepare the candidates in your input workbook: column B contains the number of elements, column C the space-separated element symbols, and column D the corresponding space-separated ratios. Ensure the elements are supported by the property and enthalpy databases. Manually transfer the generated compositions to the input workbook.

**Calculate the Model Descriptors**
Run main.m. It calls calculatefeatures.m and uses properties.xlsx, OQMD_Enthalpy.xlsx, and Orthogonalization_matrix.xlsx to calculate descriptors, including the minimum Hessian eigenvalue at the average melting temperature. The current implementation also calls SE.m.
The trained ANN requires 14 descriptors in a specific order. Ensure that the descriptor array allocated in main.m matches the number and order of the outputs from calculatefeatures.m, and that the descriptor definitions and units are consistent with those used to train the model.

**Normalize and Prepare the Inputs**
Place the 14 raw descriptor columns in Sheet1, columns I–V, then use normalization.m to prepare normalized inputs. The saved script fits min–max scaling to the current batch; new predictions should instead reuse the training normalization parameters. Save the correctly prepared inputs in Sheet2, columns I–V, for prediction. Descriptor export and transfer between sheets must be performed manually.

**Predict and Select Compositions**
Run NNprediction.m after checking its workbook name and data range. It reads the prepared inputs and calls gh_ANN.m. The output y contains one spinodal probability per alloy. Keep its row order aligned with the original compositions. To visualize screening results, plot average VEC, mixing enthalpy deviation, and minimum Hessian eigenvalue, using probability as the point color. Select promising ratios from the corresponding composition rows. A final screening plot is not generated automatically by these scripts.

**Optional Model Analysis**
SD_SS_Multi.xlsx contains 484 alloy records. Shapley_value.m trains a 14-input ANN and evaluates descriptor contributions using predict.m. These files support training and interpretation but are not required to use the exported model. SHAP analysis also requires Statistics and Machine Learning Toolbox.

**Reference**
Gong et al. (2024). Thermodynamic Descriptors for Rapid Search of Compositional Complex Spinodal Alloys with Artificial Neural Network. https://doi.org/10.1007/s44210-024-00047-x
