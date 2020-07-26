package info;

import java.util.Scanner;
import info.CheckTCP;
public class tcpLifetimeMessage {

	public static void main(String[] args) {
		
		checking();
		CheckTCP ch = new CheckTCP();
		ch.check_tcp(3, 5, 7);
	}

	private static void checking() {
		Scanner in = new Scanner(System.in);
		System.out.println("This is a prototype of the TCP Packet Lifetime massage system");
		System.out.println("Please enter The TCP Packet Lefetime (Enter a value > 0)");
		for ( int i=0; i< 5 ;i++) {
			double val = in.nextDouble();
			if (val <= 0)
				System.out.printf("\n%-5s%-15s", val, "Invalid PAcket Lifetime(Enter a value > 0)\n" );
			else if (val < 1)
				System.out.printf("\n%-5s%-15s",val, "Insuffcient Lifetime\n" );
			else if (val == 8)
				System.out.printf("\n%-5s%-15s", val, "At Least Default\n");
			else if (val > 8 && val <= 60)
				System.out.printf("\n%-5s%-15s", val, "Larger Than Default Lifetime\n");
			else if (val > 60)
				System.out.printf("\n%-5s%-15s", val, "Maximun Exceed\n");	
		}
		System.out.println("5 Packets Have Been Recorded. Message Stopped\n");
	}

}


