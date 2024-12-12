# Import the calendar module
import calendar

# User input for month and year
year = int(input("Enter year: "))
month = int(input("Enter month (1-12): "))

# Display the calendar
print(calendar.month(year, month))
