# wrapper to use callbacks.pysp_phrhosetter_callback with a different rho_cost_multiplier
import callbacks
callbacks.rho_cost_multiplier = 1.00
pysp_phrhosetter_callback = callbacks.pysp_phrhosetter_callback
