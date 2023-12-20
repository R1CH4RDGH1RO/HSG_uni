"""
Program that calculates the European/American Call/Put Option Price for a stock given user inputs.
The calculation of the European option price is based on the Black-Scholes model formula.
The calculation of the American option price is based on the Binomial model.
In the Binomial model, it is used the CRR Exact Solution approach for the calculation of u and d factors.
"""

import math

# Create Option class
class Option:
    def __init__(self, S, K, T, r, y, sigma):
        self.S = S  # Current stock price
        self.K = K  # Option strike price
        self.T = T  # Time to expiration (in years)
        self.r = r  # Risk-free interest rate
        self.y = y  # Dividend yield
        self.sigma = sigma  # Volatility of the underlying stock

# European Option class taking the Option inheritance
class EuropeanOption(Option):
    def calculate_option_price(self):
        raise NotImplementedError("Subclasses must implement the calculate_option_price method.")

# American Option class taking the Option inheritance
class AmericanOption(Option):
    def calculate_option_price(self):
        raise NotImplementedError("Subclasses must implement the calculate_option_price method.")

# European Call Option class taking the European Option inheritance
class EuropeanCallOption(EuropeanOption):
    def calculate_option_price(self):
        d1 = (math.log(self.S / self.K) + (self.r - self.y + 0.5 * self.sigma**2) * self.T) / (self.sigma * math.sqrt(self.T))
        d2 = d1 - self.sigma * math.sqrt(self.T)
        N_d1 = 0.5 * (1 + math.erf(d1 / math.sqrt(2)))
        N_d2 = 0.5 * (1 + math.erf(d2 / math.sqrt(2)))

        call_option_price = math.exp(-self.y * self.T) * self.S * N_d1 -  math.exp(-self.r * self.T) * self.K * N_d2
        return call_option_price

# European Put Option class taking the European Option inheritance
class EuropeanPutOption(EuropeanOption):
    def calculate_option_price(self):
        d1 = (math.log(self.S / self.K) + (self.r - self.y + 0.5 * self.sigma**2) * self.T) / (self.sigma * math.sqrt(self.T))
        d2 = d1 - self.sigma * math.sqrt(self.T)
        N_minus_d1 = 0.5 * (1 - math.erf(d1 / math.sqrt(2)))
        N_minus_d2 = 0.5 * (1 - math.erf(d2 / math.sqrt(2)))

        put_option_price = math.exp(-self.r * self.T) * self.K * N_minus_d2 - math.exp(-self.y * self.T) * self.S * N_minus_d1
        return put_option_price

# American Call Option class taking the American Option inheritance
class AmericanCallOption(AmericanOption):
    def calculate_option_price(self):
        num_steps = int(self.T * 252)  # One step for each trading day (252 in 1 year), taking the integer to do the for loop later
        dt = self.T / num_steps
        u = math.exp((self.r - self.y) * dt + self.sigma * math.sqrt(dt))
        d = 1/u  # CRR exact solution approach for up and down factors determination
        r = math.exp((self.r - self.y) * dt)
        p = (r - d) / (u - d)  # Risk neutral probability

        # Initialize option values at expiration
        call_values = [max(0, self.S * u**j * d**(num_steps - j) - self.K) for j in range(num_steps + 1)]

        # Iterate backward through the tree to calculate option values at each step
        for i in range(num_steps - 1, -1, -1):  # Backward approach to price the derivative (-1 step)
            call_values = [max(self.S * u**j * d**(i - j) - self.K, math.exp(-self.r * dt) * (p * call_values[j] + (1 - p) * call_values[j + 1])) for j in range(i + 1)]

        return call_values[0]

# American Put Option class taking the American Option inheritance
class AmericanPutOption(AmericanOption):
    def calculate_option_price(self):
        num_steps = int(self.T * 252)  # One step for each trading day (252 in 1 year), taking the integer to do the for loop later
        dt = self.T / num_steps
        u = math.exp((self.r - self.y) * dt + self.sigma * math.sqrt(dt))
        d = 1/u  # CRR exact solution approach for up and down factors determination
        r = math.exp((self.r -self.y) * dt)
        p = (r - d) / (u - d)  # Risk neutral probability

        # Initialize option values at expiration
        put_values = [max(0, self.K - self.S * u**j * d**(num_steps - j)) for j in range(num_steps + 1)]

        # Iterate backward through the tree to calculate option values at each step
        for i in range(num_steps - 1, -1, -1):  # Backwards approach to price the derivative (-1 step)
            put_values = [max(self.K - self.S * u**j * d**(i - j), math.exp(-self.r * dt) * (p * put_values[j] + (1 - p) * put_values[j + 1])) for j in range(i + 1)]

        return put_values[0]

# User input function to ask for the values
def get_user_input():
    # Check for valid current stock price
    while True:
        try:
            S = float(input("Enter current stock price: "))
            if S > 0:
                break
            else:
                print("Invalid input. Please enter a positive number.")
        except ValueError:
            print("Invalid input. Please enter a (positive) number.")

    # Check for valid option strike
    while True:
        try:
            K = float(input("Enter option strike price: "))
            if K > 0:
                break
            else:
                print("Invalid input. Please enter a positive number.")
        except ValueError:
            print("Invalid input. Please enter a (positive) number.")

    # Check for valid time to expiration
    while True:
        try:
            T = float(input("Enter time to expiration (in years): "))
            if T > 0:
                break
            else:
                print("Invalid input. Please enter a positive number.")
        except ValueError:
            print("Invalid input. Please enter a (positive) number.")
    
    # Check for valid risk-free interest rate
    try:
        r = float(input("Enter risk-free interest rate (in this way: 0.03 to indicate 3%), it can also be negative: "))
    except ValueError:
        print("Invalid input. Please enter a number.")

    # Check for valid dividend yield
    while True:
        try:
            y = float(input("Enter dividend yield rate (in this way: 0.03 to indicate 3%): "))
            if y > 0:
                break
            else:
                print("Invalid input. Please enter a positive number.")
        except ValueError:
            print("Invalid input. Please enter a (positive) number.")
    
    # Check for valid volatility
    while True:
        try:
            sigma = float(input("Enter volatility of the underlying stock (in this way: 0.03 to indicate 3%): "))
            if sigma > 0:
                break
            else:
                print("Invalid input. Please enter a positive number.")
        except ValueError:
            print("Invalid input. Please enter a (positive) number.")
    
    # Check for valid input for option_type
    while True:
        option_type = input("Enter the option type, E = European, A = American (it is not case sensitive): ").lower()
        if option_type in ['e', 'a']:
            break
        else:
            print("Invalid input. Please enter 'E' for European or 'A' for American.")

    # Check for valid input for call_put
    while True:
        call_put = input("Enter the option type, C = Call, P = Put (it is not case sensitive): ").lower()
        if call_put in ['c', 'p']:
            break
        else:
            print("Invalid input. Please enter 'C' for Call or 'P' for Put.")

    return S, K, T, r, y, sigma, option_type, call_put

def main():
    print("European or American Option pricing")

    # Get user input
    S, K, T, r, y, sigma, option_type, call_put = get_user_input()

    # Loop that checks the users option_type and call_put and then calculates and displays the option price
    if option_type == 'e' and call_put == 'c':
        european_call = EuropeanCallOption(S, K, T, r, y, sigma)
        european_call_price = european_call.calculate_option_price()
        print(f"European Call Option Price: {european_call_price:.4f}")

    elif option_type == 'e' and call_put == 'p':
        european_put = EuropeanPutOption(S, K, T, r, y, sigma)
        european_put_price = european_put.calculate_option_price()
        print(f"European Put Option Price: {european_put_price:.4f}")

    elif option_type == 'a' and call_put == 'c':
        american_call = AmericanCallOption(S, K, T, r, y, sigma)
        american_call_price = american_call.calculate_option_price()
        print(f"American Call Option Price: {american_call_price:.4f}")

    else:  # option_type == 'a' and call_put == 'p':
        american_put = AmericanPutOption(S, K, T, r, y, sigma)
        american_put_price = american_put.calculate_option_price()
        print(f"American Put Option Price: {american_put_price:.4f}")

# Construct that makes sure that the code will only run when the script is executed directly,
# and not when it's imported as a module into another script
if __name__ == "__main__":
    main()
