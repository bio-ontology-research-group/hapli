# In src/main.py -> configure_logging
# ...
# Remove existing handlers attached to the root logger
for handler in root_logger.handlers[:]: # <-- Removes ALL handlers from root
    root_logger.removeHandler(handler)
    handler.close()
# ...
# Add the handler to the root logger
root_logger.addHandler(handler) # Adds the tool's handler
